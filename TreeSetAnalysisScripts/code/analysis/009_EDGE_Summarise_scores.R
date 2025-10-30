# Script: 009_EDGE_Summarise_scores.R
# Description: This script summarizes EDGE (Evolutionarily Distinct and Globally Endangered) 
# scores for languages across multiple phylogenetic trees, computes rankings, 
# and merges results with metadata.
#------------------------------------------------------------------------------
# Input data:
# - EDGE_scores/ [Directory: EDGE score files, pattern "_EDGESCORES2.csv"]
# - edscores/ [Directory: Additional EDGE score files, pattern ".csv"]
# - input_data/languages_and_dialects_geoFINALupdates6.csv [CSV: Language metadata]

# Output data:
# - outputs/EDGEscores/ [Directory: Summary EDGE results]
#   - ALL_EDGE.csv [CSV: All EDGE scores]
#   - ALL_ED_wTREES.csv [CSV: EDGE scores with trees]
#   - NEW_EDGE_top_100.csv [CSV: Top 100 for 226 langs]
#   - ALL_EDGE_top_100_86928.csv [CSV: Top 100 summary]
#   - ALL_EDGE_top_100.csv [CSV: Top 100 summary with classification]

#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
library(ggplot2)     # For plotting
library(sp)          # Spatial data handling
library(rgdal)       # Geospatial data abstraction library
library(sf)          # Simple features for R
library(ggpubr)      # Publication ready plots
library(ggthemes)    # Extra themes for ggplot2
library(cowplot)     # Plot grid layouts
library(data.table)  # Fast data manipulation
library(here)        # Path handling

#------------------------------------------------------------------------------
# Part 1: Process EDGE scores from multiple trees
#------------------------------------------------------------------------------

# Get list of EDGE score files
fil1 <- list.files(here("EDGE_scores"), pattern="_EDGESCORES2.csv", full.names=TRUE)
fil2 <- list.files(here("EDGE_scores"), pattern="_EDGESCORES2.csv", full.names=FALSE)
fil3 <- read.table(text=fil2, sep="_")

# Process each file
for (j in 1:length(fil1)) {
  # Read the current file
  EDGE1 <- read.csv(fil1[j], stringsAsFactors = FALSE)
  
  # Prepare data frame
  EDGE1b <- EDGE1
  EDGE1b$tree <- fil3$V2[j]
  
  # Calculate EDGE rank and binary indicators for different rank thresholds
  EDGE1b$EDGErank <- rank(1/EDGE1b$EDGE)
  
  # Create binary indicators for ranks (1 if rank is <= threshold, 0 otherwise)
  EDGE1b$EDGErank1 <- 0; EDGE1b$EDGErank1[EDGE1b$EDGErank <= 1] <- 1
  EDGE1b$EDGErank5 <- 0; EDGE1b$EDGErank5[EDGE1b$EDGErank <= 5] <- 1
  EDGE1b$EDGErank10 <- 0; EDGE1b$EDGErank10[EDGE1b$EDGErank <= 10] <- 1
  EDGE1b$EDGErank15 <- 0; EDGE1b$EDGErank15[EDGE1b$EDGErank <= 15] <- 1
  EDGE1b$EDGErank20 <- 0; EDGE1b$EDGErank20[EDGE1b$EDGErank <= 20] <- 1
  EDGE1b$EDGErank25 <- 0; EDGE1b$EDGErank25[EDGE1b$EDGErank <= 25] <- 1
  EDGE1b$EDGErank30 <- 0; EDGE1b$EDGErank30[EDGE1b$EDGErank <= 30] <- 1
  EDGE1b$EDGErank35 <- 0; EDGE1b$EDGErank35[EDGE1b$EDGErank <= 35] <- 1
  EDGE1b$EDGErank40 <- 0; EDGE1b$EDGErank40[EDGE1b$EDGErank <= 40] <- 1
  EDGE1b$EDGErank45 <- 0; EDGE1b$EDGErank45[EDGE1b$EDGErank <= 45] <- 1
  EDGE1b$EDGErank50 <- 0; EDGE1b$EDGErank50[EDGE1b$EDGErank <= 50] <- 1
  EDGE1b$EDGErank55 <- 0; EDGE1b$EDGErank55[EDGE1b$EDGErank <= 55] <- 1
  EDGE1b$EDGErank60 <- 0; EDGE1b$EDGErank60[EDGE1b$EDGErank <= 60] <- 1
  EDGE1b$EDGErank65 <- 0; EDGE1b$EDGErank65[EDGE1b$EDGErank <= 65] <- 1
  EDGE1b$EDGErank70 <- 0; EDGE1b$EDGErank70[EDGE1b$EDGErank <= 70] <- 1
  EDGE1b$EDGErank75 <- 0; EDGE1b$EDGErank75[EDGE1b$EDGErank <= 75] <- 1
  EDGE1b$EDGErank80 <- 0; EDGE1b$EDGErank80[EDGE1b$EDGErank <= 80] <- 1
  EDGE1b$EDGErank85 <- 0; EDGE1b$EDGErank85[EDGE1b$EDGErank <= 85] <- 1
  EDGE1b$EDGErank90 <- 0; EDGE1b$EDGErank90[EDGE1b$EDGErank <= 90] <- 1
  EDGE1b$EDGErank95 <- 0; EDGE1b$EDGErank95[EDGE1b$EDGErank <= 95] <- 1
  EDGE1b$EDGErank100 <- 0; EDGE1b$EDGErank100[EDGE1b$EDGErank <= 100] <- 1
  
  # Combine results
  if(j == 1) { 
    EDGE2 <- EDGE1b 
  } else {
    EDGE2 <- rbind(EDGE2, EDGE1b)
  }
  
  print(paste("Processing file", j, "of", length(fil1)))
}

# Save combined results
# fwrite(EDGE2, file=here("outputs", "EDGEscores", "ALL_EDGE.csv"))

# Read saved data
EDGE2 <- read.csv(here("outputs", "EDGEscores", "ALL_EDGE.csv"), stringsAsFactors = FALSE)

#------------------------------------------------------------------------------
# Part 2: Summarize EDGE scores by language
#------------------------------------------------------------------------------

# Calculate mean values by species
tt <- aggregate(EDGE2[, c(5, 6, 8:(ncol(EDGE2)))], by=list(EDGE2$Species), mean, na.rm=TRUE)

# Calculate standard deviation by species
tt2 <- aggregate(EDGE2[, c(5, 6, 8:(ncol(EDGE2)))], by=list(EDGE2$Species), sd, na.rm=TRUE)

# Calculate median values by species
tt3 <- aggregate(EDGE2[, c(5, 6, 8:(ncol(EDGE2)))], by=list(EDGE2$Species), median, na.rm=TRUE)

# Calculate mean ranks for different thresholds
tt$meanEDGErank25 <- rowMeans(tt[, 5:9])
tt$meanEDGErank50 <- rowMeans(tt[, 5:14])
tt$meanEDGErank100 <- rowMeans(tt[, 5:ncol(tt)])

# Sort by meanEDGErank100 in descending order
tt <- tt[order(tt$meanEDGErank100, decreasing=TRUE), ]

# Fix specific naming issue
tt$Group.1[tt$Group.1 == "iron1242"] <- "osse1243"

#------------------------------------------------------------------------------
# Part 3: Process ED (Evolutionary Distinctiveness) scores
#------------------------------------------------------------------------------

# Get list of ED score files
fil1 <- list.files(here("edscores_equalsplits"), pattern=".csv", full.names=TRUE)

# Process each file
for (j in 1:length(fil1)) {
  # Extract tree number from filename
  tree_nr = gsub(".csv", "", basename(fil1[j]))
  tree_nr = gsub("edscores_", "", tree_nr)
  
  # Read the current file
  ED1 <- read.csv(fil1[j], stringsAsFactors = FALSE)
  ED1$tree <- tree_nr
  
  # Combine results
  if(j == 1) { 
    ED2 <- ED1 
  } else {
    ED2 <- rbind(ED2, ED1)
  }
  
  print(paste("Processing ED file", j, "of", length(fil1)))
}

# Save combined results
fwrite(ED2, file=here("outputs", "EDGEscores", "ALL_ED_wTREES.csv"))

# Read EDGE data again (this seems redundant in the original script)
ED2 <- read.csv(here("outputs", "EDGEscores", "ALL_EDGE.csv"), stringsAsFactors = FALSE)
#fwrite(ED2, file=here("outputs", "EDGEscores", "ALL_EDGE_eq.csv"))
ED2$Species[ED2$Species == "iron1242"] <- "osse1243"
#------------------------------------------------------------------------------
# Part 4: Merge EDGE scores with language metadata
#------------------------------------------------------------------------------
# I want to copy ED values from previous ED2 into ED2 where $species and $tree match
# Merge ED scores with EDGE scores
all = read.csv(here('outputs','EDGEscores','ALL_ED_wTREES.csv'))
all$glottocode[all$glottocode == "iron1242"] <- "osse1243"

# Join only on the columns needed and update the ED column in ED2
ED2 <- ED2 %>%
  select(-ED) %>%  # Remove original ED to avoid confusion
  left_join(all %>% select(glottocode, ED, tree),
            by = c("Species" = "glottocode", "tree" = "tree")) %>% relocate(ED, .after = pext)

# Calculate mean values by species
tt <- aggregate(ED2[, c(5, 6, 8:(ncol(ED2)))], by=list(ED2$Species), mean, na.rm=TRUE)
# Calculate standard deviation by species
tt2 <- aggregate(ED2[, c(5, 6, 8:(ncol(ED2)))], by=list(ED2$Species), sd, na.rm=TRUE)
# Calculate median values by species
tt3 <- aggregate(ED2[, c(5, 6, 8:(ncol(ED2)))], by=list(ED2$Species), median, na.rm=TRUE)

# Fix specific naming issue
tt$Group.1[tt$Group.1 == "iron1242"] <- "osse1243"

# Read language data
res5 <- read.csv(file=here("input_data", "languages_and_dialects_geoFINALupdates6.csv"))

# Merge EDGE scores with language names
all_names <- merge(tt, res5[, c("glottocode", "name", "iso_final")], by.x="Group.1", by.y="glottocode")

# Filter for top 100 languages
all_names2 <- all_names[all_names$EDGErank100 > 0, ]

# Save top languages
fwrite(all_names2, file=here("outputs", "EDGEscores", "NEW_EDGE_top_100.csv"))
fwrite(all_names, file=here("outputs", "EDGEscores", "ALL_EDGE_top_100.csv"))

# Convert to data.table for faster operations
setDT(all_names2)
setkey(all_names2, iso_final)
setDT(res5)
setkey(res5, iso_final)

# Clean data - filter out records with missing ISO codes
res5 <- res5[!is.na(ISO_639), ]
setkey(res5, glottocode)

# Final merge and output
setDT(tt)
setkey(tt, Group.1)
tt2b <- res5[tt]

# Write final results
fwrite(tt2b, file=(here("outputs", "ALL_EDGE_top_100_86928.csv")))