# Script: 012_LTT_data.R
# Description: Creates data for Lineage Through Time (LTT) plots for Figure 1.
#              Processes phylogenetic trees by continent and generates LTT metrics.

#------------------------------------------------------------------------------
# Input data:
# - input_data/languages_and_dialects_geoFINALupdates6.csv [CSV: Language metadata]
# - all_trees/ [Directory: Phylogenetic trees by continent]

# Output data:
# - outputs/Fig_1/LTT_Plot1_Data.csv [CSV: LTT plot data]

#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
library(treespace)    # Tree space analysis
library(phangorn)     # Phylogenetic analysis
library(phytools)     # Phylogenetic tools
library(phyclust)     # Phylogenetic clustering
library(ggsci)        # Scientific journal color palettes
library(ggtree)       # Visualization of phylogenetic trees
library(paleotree)    # Paleontological and phylogenetic tools
library(diversitree)  # Comparative phylogenetic methods
library(ape)          # Analyses of Phylogenetics and Evolution
library(data.table)   # Enhanced data frame
library(rgdal)        # Geospatial data abstraction
library(sp)           # Classes and methods for spatial data
library(here)         # Path handling

#------------------------------------------------------------------------------
# Part 1: Load language data
#------------------------------------------------------------------------------
# Read language data with dialect information
res5 <- read.csv(file = here("input_data", 
                            "languages_and_dialects_geoFINALupdates6.csv"), 
                stringsAsFactors = FALSE)

#------------------------------------------------------------------------------
# Part 2: Generate LTT data by continent
#------------------------------------------------------------------------------
# Define continents for analysis
all_conts <- c("Global", "Oceania", "Africa", "Eurasia",
               "North-America", "South-America")

# Define sample size (using 1000 trees per continent)
samp1 <- 1:1000

# Initialize data storage
trees <- NULL

# Loop through each continent
for (kk in 1:length(all_conts)) {
  # Get trees for current continent
  tt <- list.files(here("all_trees"), pattern = all_conts[kk], full.names = TRUE)
  
  # Loop through tree samples
  for(ii in samp1) {
    # Read phylogenetic tree
    tr <- read.tree(tt[ii])
    
    # Check for tips not in language data
    tr$tip.label[!tr$tip.label %in% res5$id]
    
    # Create Lineage Through Time data (without plotting)
    ltt1 <- ltt(tr, plot = FALSE)
    
    # Normalize time values
    ltt1$times <- ltt1$times - max(ltt1$times)
    
    # Calculate log of lineage counts
    ltt1$ltt2 <- log(ltt1$ltt)
    
    # Store LTT metrics in data frame
    resx1 <- data.frame(
      tree_id_1_901 = ii,
      cont = all_conts[kk],
      ltt = ltt1$ltt,
      logltt = ltt1$ltt2,
      times = ltt1$times,
      gamma = ltt1$gamma,
      p = ltt1$p,
      max_height = max(nodeHeights(tr))
    )
    
    # Append to results
    if(is.null(trees) == TRUE) {
      resx2 <- resx1
      trees <- 1
    } else {
      resx2 <- rbind(resx2, resx1)
    }
  }
  
  print(paste("Processed continent:", kk, "-", all_conts[kk]))
}

#------------------------------------------------------------------------------
# Part 3: Process and save results
#------------------------------------------------------------------------------
# Store LTT data
ltt1 <- resx2

# Standardize continent names
ltt1$cont[ltt1$cont == "North-America"] <- "North America"
ltt1$cont[ltt1$cont == "South-America"] <- "South America"

# Create output directory if it doesn't exist
if(!dir.exists(here("outputs", "Fig_1"))) {
  dir.create(here("outputs", "Fig_1"), recursive = TRUE)
}

# Save LTT data
write.csv(ltt1, here("outputs", "Fig_1", "LTT_Plot1_Data.csv"), row.names = FALSE)