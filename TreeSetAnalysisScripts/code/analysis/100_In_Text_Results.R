#ARTUR - retrieval of data to input in the manuscript

# Load the necessary libraries
library(tidyverse)
library(readxl)
library(here)
library(magrittr)
library(countrycode)
library(rgdal)
library(data.table)
library(coda)
library(ape)       # For phylogenetic analysis
library(here)      # For file path handling
library(here)
library(data.table)
library(dplyr)
library(magrittr)


# Phylogenetic diversity of the world’s languages -------------------------------
## Summary of tree stats ----
grps = c("Global", "Oceania", "Africa", "Eurasia", "South America", "North America")

# Initialize an empty list to store results for both variables
results_list <- list()

# Variables to analyze
variables <- c("tree_length", "crown_age")
df <- read.csv(here("outputs","All_Trees_Summary.csv"))
# Iterate over variables and groups

for (var in variables) {
  for (grp in grps) {
    # Load the data - summary tree statistics - FULL 
    gr1 <- df %>% filter(group == grp)
    gr1_median <- median(gr1[[var]])
    
    # Load the HPDI data 
    gr1_hpd <- HPDinterval(as.mcmc(gr1[[var]]), prob = 0.95)
    gr1_95lower <- gr1_hpd[1, "lower"]
    gr1_95upper <- gr1_hpd[1, "upper"]
    
    # Format output
    lower_95 <- ifelse(var == "crown_age", signif(gr1_95lower, digits = 4), round(gr1_95lower / 1000, 2))
    upper_95 <- ifelse(var == "crown_age", signif(gr1_95upper, digits = 5), round(gr1_95upper / 1000, 2))
    median_tl <- ifelse(var == "crown_age", signif(gr1_median, digits = 4), signif(gr1_median / 1000, digits = 4))
    range_str <- paste0(lower_95, "-", upper_95)
    
    # Add output to the results list
    results_list[[paste(grp, var, sep = "_")]] <- data.frame(
      group = grp,
      variable = var,
      lower_95_HPDI = lower_95,
      upper_95_HPDI = upper_95,
      median_value = median_tl,
      range = range_str
    )
  }
}

# Combine all results into a single data frame
results_df <- do.call(rbind, results_list)

# View the results
print(results_df)

# Save the results to a CSV file if needed
write.csv(results_df, here("outputs","in_text_data", "summary_tree_stats.csv"), row.names = FALSE)

# Summarise tip ages  ----------------------------------------------------- 
# Define dataset source ("old" or "new")
tree_set <- "new"  # Change to "new" if needed

# Set folder and output filename based on tree_set
if (tree_set == "old") {
  folder_path <- here("input_data", "all_trees")
  output_file <- "summary_ages_old_trees.csv"
} else if (tree_set == "new") {
  folder_path <- here("all_trees_global")
  output_file <- "summary_ages_new_trees.csv"
} else {
  stop("Invalid tree_set value. Use 'old' or 'new'.")
}

# Get tree files from the specified folder
tree_files <- list.files(folder_path, pattern = "\\.tree$", full.names = TRUE)

# Report the number of trees found
cat(length(tree_files), "trees found.\n")

# Initialize vectors for global statistics and per-tree metrics
all_tips_all_trees <- numeric(0)  # Stores all tip lengths across all trees
all_tree_medians   <- numeric(length(tree_files))  # Median per tree
percent_in_5k      <- numeric(length(tree_files))  # % in last 5000 years per tree

# Process each tree
for (i in seq_along(tree_files)) {
  
  # Read and format tree
  tree2 <- read.tree(tree_files[i])
  tree  <- as.phylo(tree2)
  
  # Extract tip labels and edge lengths
  tips   <- tree$tip.label
  edges  <- tree$edge.length
  
  # Get number of tips and identify tip rows
  n_tips   <- length(tips)
  tip_rows <- which(tree$edge[, 2] <= n_tips)
  
  # Calculate tip lengths (ages)
  tip_lengths <- edges[tip_rows]
  
  # Append these tip lengths to the global vector
  all_tips_all_trees <- c(all_tips_all_trees, tip_lengths)
  
  # (Optional) For each tree, calculate median tip length
  #all_tree_medians[i] <- median(edges, na.rm = TRUE)
  all_tree_medians[i] <- median(tip_lengths, na.rm = TRUE)
  
  # Percentage of tips (languages) within the last 5000 years for this tree
  percent_in_5k[i] <- sum(tip_lengths*1000 <= 5000, na.rm = TRUE) / length(tip_lengths) * 100
}

## 1. Age below which 95% of languages fall (global statistic) ------
cutoff_95_all <- quantile(all_tips_all_trees, probs = 0.95, na.rm = TRUE)
print(cutoff_95_all*1000)

median(all_tips_all_trees)*1000
mean((all_tips_all_trees)*1000)

## 2. Compute 95% HPD interval for the % of languages in the last 5000 years ----
percent_in_5k_mcmc <- as.mcmc(percent_in_5k)
percent_in_5k_hpd  <- HPDinterval(percent_in_5k_mcmc, prob = 0.95)

# (Optional) Summaries of the per-tree median
medians_mcmc <- as.mcmc(all_tree_medians)
medians_hpd  <- HPDinterval(medians_mcmc, prob = 0.95)
print(medians_hpd*1000)

# (Optional) Summaries of the per-tree median
all_tips_mcmc <- as.mcmc(all_tips_all_trees)
all_tips_hpd  <- HPDinterval(all_tips_mcmc, prob = 0.95)
print(all_tips_hpd*1000)


## Combine results into a summary data frame -----
summary_lang_ages <- data.frame(
  mean_tree_median  = mean(all_tree_medians, na.rm = TRUE)*1000,
  hpd_lower_95_median = medians_hpd[1]*1000,
  hpd_upper_95_median = medians_hpd[2]*1000,
  cutoff_95_all      = cutoff_95_all*1000,
  mean_percent_5k    = mean(percent_in_5k, na.rm = TRUE),
  hpd_lower_95_5k    = percent_in_5k_hpd[1],
  hpd_upper_95_5k    = percent_in_5k_hpd[2]
)

cat("Summary of language ages:\n")
cat("Mean tree median: ", round(summary_lang_ages$mean_tree_median, 1), " years\n")
cat("95% HPD interval for tree median: ", 
    round(summary_lang_ages$hpd_lower_95_median, 1), " - ", 
    round(summary_lang_ages$hpd_upper_95_median, 1), " years\n")
cat("Cutoff below which 95% of languages fall: ", 
    round(summary_lang_ages$cutoff_95_all, 1), " years\n")

# Print results
print(round(summary_lang_ages,1))

# Save summary to CSV
write.csv(
  summary_lang_ages,
  file = here("outputs","in_text_data", output_file),
  row.names = FALSE
)

# ========================================================================== #

# Description:
# This R script processes and analyzes summary statistics for growth rates (median, minimum, and maximum) across tree species.
# It consolidates data from multiple CSV files stored in the 'summary_ARTUR' directory, merging them into a single dataset.
# The script computes key descriptive statistics such as median and mean values for each growth rate metric.
# Additionally, it calculates 95% Highest Density Posterior Intervals (HDPI) for each metric using the 'coda' package.

# Key Steps:
# 1. Import necessary libraries and set up paths to input and output files.
# 2. Load and combine data from multiple summary files into a unified data table.
# 3. Compute median and mean values for growth rates across datasets.
# 4. Use Bayesian analysis with the 'coda' package to estimate 95% HDPI intervals.
# 5. Generate summary tables and export results to CSV files for further analysis.

# Outputs:
# 1. Combined dataset saved as 'Global_trees_summary_ARTUR.csv'.
# 2. Statistical summaries including medians and 95% HDPI intervals saved as 'DR_rates_across_tips.csv'.

# 1. List all summary files
files <- list.files(
  path = here("summary"),
  pattern = "^summary_.*\\.csv$",
  full.names = TRUE
)

# 2. Read & combine
df_all <- rbindlist(lapply(files, fread))

# 3. Inspect
dim(df_all)
df_global = df_all %>% filter(type=="all")
#write.csv(df_global, file = here("outputs", "Global_trees_summary.csv"), row.names = FALSE)

# Median/Mean Median
median(df_all$medianDR)
mean(df_all$medianDR)

## 95% HDPI for medianDR

# Median Minimum
median(df_all$minDR)
## 95% HDPI for minDR

# Median Maximum
median(df_all$maxDR)
## 95% HDPI for maxDR

# Load the 'coda' library
library(coda)

# Function to calculate 95% HDPI using 'coda'
calculate_hdpi_coda <- function(data_vector) {
  # Convert the data vector to an MCMC object
  mcmc_obj <- as.mcmc(data_vector)
  
  # Calculate the HDPI interval
  hdpi <- HPDinterval(mcmc_obj, prob = 0.95)
  
  return(hdpi)
}

# Calculate 95% HDPI for 'medianDR'
hdpi_medianDR <- calculate_hdpi_coda(df_all$medianDR)

# Calculate 95% HDPI for 'minDR'
hdpi_minDR <- calculate_hdpi_coda(df_all$minDR)

# Calculate 95% HDPI for 'maxDR'
hdpi_maxDR <- calculate_hdpi_coda(df_all$maxDR)

# Display the results
cat("95% HDPI for medianDR:\n")
print(hdpi_medianDR)
cat("\n95% HDPI for minDR:\n")
print(hdpi_minDR)
cat("\n95% HDPI for maxDR:\n")
print(hdpi_maxDR)

# Using 'coda'
summary_table_coda <- data.frame(
  Statistic = c("medianDR", "minDR", "maxDR"),
  Median = c(median(df_all$medianDR, na.rm = TRUE),
             median(df_all$minDR, na.rm = TRUE),
             median(df_all$maxDR, na.rm = TRUE)),
  HDPI_Lower = c(hdpi_medianDR[1],
                 hdpi_minDR[1],
                 hdpi_maxDR[1]),
  HDPI_Upper = c(hdpi_medianDR[2],
                 hdpi_minDR[2],
                 hdpi_maxDR[2])
)

print(summary_table_coda)
#write.csv(summary_table_coda, file = here("outputs","in_text_data","DR_rates_across_tips.csv"), row.names = FALSE)

# Family-specifci DF rates factors ----------------------------------------------

## Option one: calculate median mean DR for each family, and then calculate ratio
##  of min/max of these families
fam = df_all %>% filter(type=="top27families")

# 1. For each family & each posterior tree, compute the mean of meanDR
fam_tree_means <- fam %>%
  group_by(group) %>%
  summarise(family_meanDR = median(HmeanDR), .groups = "drop")

print(fam_tree_means %>% filter(group %in% c("Quechuan", "Uto-Aztecan")))

# 3. Compute the ratio: (maximum across families) / (minimum across families)
factor_diff <- max(fam_tree_means$family_meanDR) / 
  min(fam_tree_means$family_meanDR)

factor_diff

## Option two: calculate ratio for each family, and then calculate median of these ratios
# So here we check by what order of magnitude does DR differ within each family and 
# then what is the average ratio of these for an family

# 1. For each combination of a tree & family, compute the mean of ratio of min/max DR
fam_tree_min_max <- fam %>% mutate(ratio = maxDR/minDR)

# 2. For each family, compute the median of the ratio
fam_tree_min_max_sum <- fam_tree_min_max %>%
  group_by(group) %>%
  summarise(mean_ratio = median(ratio), .groups = "drop")

# 3. Compute the median of the mean ratios across families
median(fam_tree_min_max_sum$mean_ratio)


