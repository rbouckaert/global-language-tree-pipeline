# Load the necessary libraries
library(tidyverse)
library(readxl)
library(here)
library(magrittr)
library(countrycode)
library(rgdal)
library(data.table)
library(coda)
library(sp)
library(dplyr)
library(ggplot2)
library(maps)        # For map_data("world")
library(viridis)  
library(psych)
library(sf)
library(rnaturalearth)

#------------------------------------------------------------------------------#
#             Supplementary Table S1   - country codes                  ------
#------------------------------------------------------------------------------#


# Supplementary Table 1 with PD/country and #of languages per country
## Get data  --------------------------------------/

# Load ED scores
EDGE2 <- read.csv(here("outputs","EDGEscores", "ALL_ED_wTREES.csv"))

tmp = EDGE2 %>% group_by(glottocode) %>% summarise(ED = mean(ED, na.rm = TRUE))

# Load country codes
res5 <- read.csv(file = here("input_data", "languages_and_dialects_geoFINALupdates6.csv"))
res5$glottocode[res5$glottocode=="osse1243"] <- "iron1242"

## Incorporate tree-to-tree uncertainty ----- ##

# Here we work with the individual ED values per tree.
# It is assumed that EDGE2 contains a column "tree" identifying each tree.
# Merge EDGE2 with res5 to assign country codes to each language for each tree.
all_names_tree <- merge(EDGE2, res5[, c("glottocode", "name", "iso_final", "country_ids")],
                        by = "glottocode")

# Expand the 'country_ids' column so that each language associated with multiple countries is properly represented.
all_names_tree_expanded <- all_names_tree %>%
  mutate(country_list = strsplit(as.character(country_ids), " ")) %>%
  unnest(country_list)

# For each tree, sum the ED scores of all languages in each country.
ed_country_tree <- all_names_tree_expanded %>%
  group_by(country_list, tree) %>%
  summarise(total_ED_tree = sum(ED, na.rm = TRUE), .groups = 'drop')

# calculate:
# - mean_total_ED: The average total ED across trees.
# - sd_total_ED: The standard deviation of the total ED across trees, representing tree-to-tree uncertainty.
ed_country_uncertainty <- ed_country_tree %>%
  group_by(country_list) %>%
  summarise(mean_total_ED = mean(total_ED_tree, na.rm = TRUE),
            sd_total_ED = sd(total_ED_tree, na.rm = TRUE),
            .groups = 'drop')

# Add full country names
ed_country_uncertainty$country_name <- countrycode(ed_country_uncertainty$country_list,
                                                   origin = "iso2c",
                                                   destination = "country.name")

# Order the results by the mean total ED and view the table
ed_country_uncertainty <- ed_country_uncertainty %>% arrange(desc(mean_total_ED))
print(ed_country_uncertainty)

ed_country_uncertainty <- ed_country_uncertainty %>%
  left_join(langs_per_country, by = "country_list") %>%
  mutate(mean_ED_per_language = mean_total_ED / num_languages) %>% 
  select(-langs)

# Save the results to a CSV file if needed
write.csv(ed_country_uncertainty, here("outputs","Table_S1", "summed_EqSplits_ED_per_country_TabS1_v2.csv"), row.names = FALSE)


#------------------------------------------------------------------------------#
#  Supplementary Table S1 - Method with centroids and tree-to-tree uncertainty -----
#------------------------------------------------------------------------------#
library(tidyverse)
library(here)
library(sf)
library(countrycode)
library(rnaturalearth)  # For ne_countries()

# Load tree-specific ED scores
EDGE2 <- read.csv(here("outputs","EDGEscores", "ALL_ED_wTREES.csv"))

# Correct glottocode if needed
EDGE2$glottocode[EDGE2$glottocode == "osse1243"] <- "iron1242"

# Load language spatial data (centroids)
langa <- read.csv(here("input_data", "languages_and_dialects_geoFINALupdates6.csv"),
                  stringsAsFactors = FALSE)
langa <- langa[!is.na(langa$latitude) & !is.na(langa$longitude), ]
langa$glottocode[langa$glottocode == "osse1243"] <- "iron1242"

# Convert langa to an sf object (each row represents a language centroid)
langa_sf <- st_as_sf(langa, coords = c("longitude", "latitude"), crs = 4326)

# Load country boundaries as an sf object
countries <- ne_countries(scale = "medium", returnclass = "sf")

# Transform langa to match country CRS and perform the spatial join (once)
langa_sf <- st_transform(langa_sf, st_crs(countries))
langa_joined <- st_join(langa_sf, countries, join = st_within)

# Remove languages that did not fall within any country
langa_joined <- langa_joined %>% filter(!is.na(admin))

# Reduce to a unique language-to-country table
language_country <- langa_joined %>%
  st_set_geometry(NULL) %>%  # drop geometry to save memory
  select(glottocode, country = admin) %>%
  distinct()

# Calculate the number of languages per country
language_count <- language_country %>%
  group_by(country) %>%
  summarise(num_languages = n(), .groups = "drop")

# Merge the country assignment with tree-specific ED scores.
EDGE2_country <- merge(EDGE2, language_country, by = "glottocode")

# For each country and tree, sum the ED values.
ed_country_tree <- EDGE2_country %>%
  group_by(country, tree) %>%   # "tree" identifies the tree-specific ED value
  summarise(total_ED_tree = sum(ED, na.rm = TRUE), .groups = "drop")

# For each country, calculate the mean and SD of total ED across trees.
ed_country_uncertainty <- ed_country_tree %>%
  group_by(country) %>%
  summarise(mean_total_ED = mean(total_ED_tree, na.rm = TRUE),
            sd_total_ED = sd(total_ED_tree, na.rm = TRUE),
            .groups = "drop")

# Merge language count into the summary
ed_country_uncertainty <- ed_country_uncertainty %>%
  left_join(language_count, by = "country")

# Calculate mean ED per language for each country
ed_country_uncertainty <- ed_country_uncertainty %>%
  mutate(mean_ED_per_language = mean_total_ED / num_languages)

# Optionally, add full country names using countrycode:
ed_country_uncertainty$country_name <- countrycode(ed_country_uncertainty$country,
                                                   origin = "iso2c",
                                                   destination = "country.name")

# Order and view the results
ed_country_uncertainty <- ed_country_uncertainty %>% arrange(desc(mean_total_ED))
print(ed_country_uncertainty)

# Save the results if desired
write.csv(ed_country_uncertainty,
          here("outputs","Table_S1","TableS1_ED_per_country_centroids_SD.csv"),
          row.names = FALSE)
