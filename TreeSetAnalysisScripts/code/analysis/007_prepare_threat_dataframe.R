#-------------------------------------------------------------------------------#
#                     Threat Data Preparation Script                            #
#-------------------------------------------------------------------------------#
# Purpose:                                                                      
#   This script prepares a comprehensive dataframe containing language threat   
#   data from multiple sources (EGIDS and Glottolog). It combines geographical  
#   information, language status, and spatial attributes to create a processed  
#   threat dataset for further analysis.                                        

#-------------------------------------------------------------------------------

# Input data:
# - input_data/languages_and_dialects_geoFINALupdates6.csv [CSV: Language metadata]
# - input_data/islands_transformed.r [RData: Island polygons]
# - input_data/threat_data1.csv [CSV: Threat assessment data]
# - input_data/langa/langa.shp [Shapefile: Language polygons]

# Output data:
# - input_data/processed_threat_data_frame.csv [CSV: Processed threat data]

#-------------------------------------------------------------------------------
#                       Load Required Libraries                                 
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
library(rworldmap)  # For world map data

#-------------------------------------------------------------------------------#
#                      Load and Prepare Base Data                               #
#-------------------------------------------------------------------------------#
# Load world map data for spatial analysis
wrld_simpl2 <- getMap()

# Load main language data file
res5 <- fread(file = here("input_data", "languages_and_dialects_geoFINALupdates6.csv"))
# Adjust subsistence index (0-based to 1-based conversion)
res5$`DPLACE subsistence` <- res5$`DPLACE subsistence` - 1

#-------------------------------------------------------------------------------#
#                      Add Geographical Information                             #
#-------------------------------------------------------------------------------#
# Convert language data to spatial points
rr <- res5
coordinates(rr) <- ~longitude + latitude
projection(rr) <- projection(wrld_simpl2)

# Extract region information from world map
rr2 <- over(rr, wrld_simpl2)

# Add region information to main dataset
res5 <- cbind(res5, rr2)
res5[, "REGION" := GEO3major]
res5[, "SUBREGION" := IMAGE24]

# Load island data and determine which languages are on islands
load(file = here("input_data", "islands_transformed.r"))
bi2$dummy <- 1
rr3 <- over(rr, bi2)
rr3$dummy[is.na(rr3$dummy)] <- 0

# Add island status to main dataset
res5[, island := rr3$dummy]
setkey(res5, iso_final)

#-------------------------------------------------------------------------------#
#                        Merge Threat Data                                      #
#-------------------------------------------------------------------------------#
# Load threat assessment data
lang_threat <- fread(here("input_data", "threat_data1.csv"))
setkey(lang_threat, "ISO")

# Filter to 40-year projection
lang_threat2 <- lang_threat[Years == 40, ]

# Rename columns for clarity
names(lang_threat2)[2:9] <- paste0("X", names(lang_threat2))[2:9]

# Merge threat data with main dataset
res5 <- lang_threat2[res5]
setkey(res5, ISO)

#-------------------------------------------------------------------------------#
#                     Add Language Area Information                             #
#-------------------------------------------------------------------------------#
# Load language shape files
langa <- readOGR(here("input_data", "langa", "langa.shp"), "langa")

# Extract shape data
langa_data <- setDT(langa@data)

# Select columns related to language area
cols_chosen <- c("Shp_Lng", "Shap_Ar")

# Calculate mean area by language
agg1 <- langa_data[order(LANG_IS), lapply(.SD, mean), by = .(LANG_IS), 
                  .SDcols = cols_chosen]
setkey(agg1, LANG_IS)

# Merge area data with main dataset
res5 <- agg1[res5]
res5[, c("V21", "V22") := NULL]

#-------------------------------------------------------------------------------#
#                     Calculate Threat Indices                                  #
#-------------------------------------------------------------------------------#
# Convert EGIDS threat levels to numeric values
res5$threatn <- res5$EGIDS
res5$threatn <- gsub("a", "", res5$threatn)
res5$threatn <- gsub("b", "", res5$threatn)
res5$threatn <- gsub("x", "", res5$threatn)
res5$threat_EGIDS <- as.numeric(res5$threatn)

# Convert Glottolog status to numeric threat levels
res5$threatn2 <- res5$status
res5$threatn2 <- gsub("critically endangered", "5", res5$threatn2)
res5$threatn2 <- gsub("definitely endangered", "3", res5$threatn2)
res5$threatn2 <- gsub("severely endangered", "4", res5$threatn2)
res5$threatn2 <- gsub("vulnerable", "2", res5$threatn2)
res5$threatn2 <- gsub("extinct", "6", res5$threatn2)
res5$threatn2 <- gsub("safe", "1", res5$threatn2)
res5$threat_glot <- as.numeric(res5$threatn2)

#-------------------------------------------------------------------------------#
#                     Export Processed Dataset                                  #
#-------------------------------------------------------------------------------#
# Write the final processed dataframe to CSV
write.csv(res5, file = here("input_data", "processed_threat_data_frame.csv"))
