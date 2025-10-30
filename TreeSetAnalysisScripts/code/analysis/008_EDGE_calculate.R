#-------------------------------------------------------------------------------#
#                    EDGE Score Calculation Script                              #
#-------------------------------------------------------------------------------#
# Purpose:                                                                      #
#   This script calculates Evolutionary Distinctive and Globally Endangered     
#   (EDGE) scores for languages based on phylogenetic trees and threat data.   
#   It processes multiple trees and creates EDGE scores for each, also          
#   calculating potential phylogenetic diversity loss under extinction          
#   scenarios.                                                                 
#-------------------------------------------------------------------------------
#                       Input and Output Data Summary                          

# Input data:
# - all_trees/ [Directory: Phylogenetic trees, files with pattern "Global"]
# - input_data/processed_threat_data_frame.csv [CSV: Processed threat data]
# - code/EDGE2.R [R script: EDGE score calculation function]

# Output data:
# - EDGE_scores/ [Directory: EDGE scores and PD loss results]
#   - [tree_name]_EDGESCORES2.csv [CSV: EDGE scores]
#   - [tree_name]_PDLOST2.csv [CSV: PD loss calculations]

#-------------------------------------------------------------------------------#
#                       Load Required Libraries                                 #
#-------------------------------------------------------------------------------#
library(ape)        # For phylogenetic analysis
library(caper)      # For comparative analysis
library(rgdal)      # For spatial data handling
library(sp)         # For spatial data structures
library(data.table) # For data manipulation
library(raster)     # For raster data handling
library(mice)       # For missing data imputation
library(phylobase)  # For phylogenetic data structures
library(here)       # For file path management

#-------------------------------------------------------------------------------#
#                       Source External Functions                               #
#-------------------------------------------------------------------------------#
# EDGE2 calculation function written by Rikki Gumbs
source(here("code", "analysis", "EDGE2.R"))

#-------------------------------------------------------------------------------#
#                       Load and Prepare Input Data                             #
#-------------------------------------------------------------------------------#
# Read in all phylogenetic trees
tree_files <- list.files(here("all_trees"), pattern="Global", full.names=TRUE)

# Load prepared language threat data
threat_data <- fread(file=here("input_data", "processed_threat_data_frame.csv"))

# Subset needed columns for EDGE calculation
# - glottocode: Language identifier
# - threat_glot: Threat level from Glottolog
# - threat_EGIDS: Threat level from EGIDS
# - L1_Users: Number of L1 speakers (log-transformed)
# - Shap_Ar: Language geographical area (log-transformed)
# - X10: Extinction probability
selected_data <- threat_data[, c("glottocode", "threat_glot", "threat_EGIDS", 
                               "L1_Users", "Shap_Ar", "X10")]

# Apply log transformations to numerical variables
selected_data$L1_Users <- log(selected_data$L1_Users + 1)
selected_data$Shap_Ar <- log(selected_data$Shap_Ar + 1)

#-------------------------------------------------------------------------------#
#                       Impute Missing Values                                   #
#-------------------------------------------------------------------------------#
# Use multiple imputation (CART method) to handle missing data
# Creating 5 imputed datasets with 50 iterations
imputation_model <- suppressWarnings(mice(selected_data, m=5, maxit=50, 
                                        method="cart", print=FALSE))

# Extract the 5 complete datasets
imputed_data_1 <- complete(imputation_model, 1)
imputed_data_2 <- complete(imputation_model, 2)
imputed_data_3 <- complete(imputation_model, 3)
imputed_data_4 <- complete(imputation_model, 4)
imputed_data_5 <- complete(imputation_model, 5)

# Calculate mean values across all imputed datasets
final_data <- data.frame(
  glottocode = imputed_data_1$glottocode,
  (imputed_data_1[, 2:6] + imputed_data_2[, 2:6] + 
   imputed_data_3[, 2:6] + imputed_data_4[, 2:6] + 
   imputed_data_5[, 2:6]) / 5
)

# Round EGIDS threat values to integers
final_data$threat_EGIDS <- round(final_data$threat_EGIDS, 0)

# Recode EGIDS level 10 to 9 (merging highest threat categories)
final_data$threat_EGIDS[final_data$threat_EGIDS == 10] <- 9

#-------------------------------------------------------------------------------#
#                       Prepare Extinction Probability Data                     #
#-------------------------------------------------------------------------------#
# Rename extinction probability column for clarity
names(final_data)[6] <- "ext_prob"

# Convert to standard data frame
extinction_data <- as.data.frame(final_data)

# Ensure minimum extinction probability for calculation
min_probability <- min(extinction_data$ext_prob[extinction_data$ext_prob != 0])
extinction_data$ext_prob[extinction_data$ext_prob <= 0] <- min_probability

# Prepare final dataset with species ID and extinction probability
edge_input <- extinction_data[, c("glottocode", "ext_prob")]
names(edge_input) <- c("species", "pext")

# Fix specific language code issue
edge_input$species[edge_input$species == "osse1243"] <- "iron1242"

# Randomize tree processing order
tree_files <- sample(tree_files)

#-------------------------------------------------------------------------------#
#                       Calculate EDGE Scores                                   #
#-------------------------------------------------------------------------------#
for (i in 1:(length(tree_files))) {
  # Generate output filenames
  edge_output_file <- gsub("all_trees", "EDGE_scores", tree_files[i], fixed=TRUE)
  edge_output_file <- gsub(".tree", "_EDGESCORES2.csv", edge_output_file, fixed=TRUE)
  
  pd_output_file <- gsub(".csv", "_PDLOST2.csv", edge_output_file, fixed=TRUE)
  
  # Skip if output file already exists
  if(edge_output_file %in% list.files(here("EDGE_scores"), full.names=TRUE)) {
    next
  }

  # Load phylogenetic tree
  tree <- read.tree(tree_files[i])
  
  # Match species in the tree to extinction data
  species_data <- edge_input[match(tree$tip.label, edge_input$species), ]
  
  # Calculate EDGE scores
  edge_results <- EDGE.2.calc(tree=tree, pext=species_data)
  
  # Extract EDGE scores table
  edge_scores <- edge_results[[1]]
  
  #---------------------------------------------------------------------------#
  #            Calculate Phylogenetic Diversity (PD) Loss                     #
  #---------------------------------------------------------------------------#
  # Calculate total phylogenetic diversity in the tree
  total_pd <- sum(tree$edge.length)
  
  # Get languages with high threat levels (EGIDS > 6)
  threatened_languages <- na.omit(final_data$glottocode[final_data$threat_EGIDS > 6])
  
  # Calculate PD loss if threatened languages go extinct
  tree_without_threatened <- drop.tip(tree, 
                                      tree$tip.label[tree$tip.label %in% threatened_languages])
  pd_loss_threatened <- total_pd - sum(tree_without_threatened$edge.length)
  
  # Calculate PD loss for a random extinction scenario (same number of languages)
  random_languages <- sample(tree$tip.label, length(threatened_languages))
  tree_random_extinction <- drop.tip(tree, 
                                    tree$tip.label[tree$tip.label %in% random_languages])
  pd_loss_random <- total_pd - sum(tree_random_extinction$edge.length)
  
  # Prepare and save PD loss results
  pd_results <- data.frame(
    tree = i,
    total_pd = total_pd,
    total_lost_threat = pd_loss_threatened,
    total_lost_random = pd_loss_random
  )
  
  # Write results to output files
  write.csv(pd_results, file=pd_output_file)
  write.csv(edge_scores, file=edge_output_file)
  
  print(paste("Processed tree", i, "of", length(tree_files)))
}
