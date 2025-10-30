# Script: 016_BISSE_summary.R
# Description: Summarises and visualises BiSSE/HiSSE model results from simulation and real data.
#------------------------------------------------------------------------------
# Input data:
# - BISSE_runs/ [Directory: BiSSE analysis results, pattern "ISSE_runs3.csv"]
# - input_data/bisse_param_names.r [RData: Parameter names for BiSSE]
#
# Output data:
# - outputs/ [Directory: Contains summary and visualisation files]
#   - bisse_results1.csv [CSV: Combined BiSSE run results]
#   - DR_test_summary.csv [CSV: DR test summary statistics]
#   - AIC_test_summary.csv [CSV: AIC summary statistics]
#   - bisse_parameter_differences.csv [CSV: Parameter differences between states]
#   - bisse_parameter_directions.csv [CSV: Proportion of runs with higher State 1 parameters]
#   - bisse_null_AIC_proportions.csv [CSV: Proportion of runs where BiSSE outperforms null model]
#------------------------------------------------------------------------------

#######################################
## SECTION 1: Load Libraries & Setup ##
#######################################

# Load required libraries for analyses and plotting
library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(ggridges)
library(data.table)
library(rethinking)
library(forcats)
library(here)  # for path management

#######################################
## SECTION 2: Read BISSE Run Files  ##
#######################################

# List all BISSE run CSV files from the BISSE_runs folder
lf1 <- list.files(here("BISSE_runs"), pattern = "ISSE_runs3.csv", full.names = TRUE)

# Extract continent information from filenames
cont1 <- gsub(paste0(here("BISSE_runs"), "/es4_tree_"), "", lf1, fixed = TRUE)
cont1 <- gsub("_BISSE_runs3.csv", "", cont1, fixed = TRUE)
cont2 <- read.table(text = cont1, sep = "_")
cont3 <- cont2$V4

# Read each BISSE run file, add continent/location, and combine into one data frame
for (i in seq_along(lf1)) {
  lf2 <- read.csv(lf1[i], stringsAsFactors = FALSE)
  lf2$location <- cont3[i]
  
  if (i == 1) {
    lf3 <- lf2
  } else {
    lf3 <- rbind(lf3, lf2)
  }
  
  print(i)
}

#######################################
## SECTION 3: Process Tree/Location ##
#######################################

# Read and parse the 'tree' column for additional information
dat1 <- read.table(text = lf3$tree, sep = "_")
names(dat1) <- c("tree2", "number", "xxxx", "type")
lf3 <- cbind(lf3, dat1)

# Write the combined results out to file
fwrite(lf3, file = here("outputs", "bisse_results1.csv"))

#######################################
## SECTION 4: Load Parameters & Test Types ##
#######################################

# Load parameter names object (assumed to be saved in 'xxx')
load(here("input_data", "bisse_param_names.r"))

# Identify test types available in the data
types <- unique(lf3$name)
print(types)

#######################################
## SECTION 5: Subset & Process DR Test ##
#######################################

# Subset for DR_test and clean variable names
DR1 <- lf3[lf3$name == "DR_test", ]
DR1$variable <- gsub("Marriage_Endogamous", "Endogamous Marriage", DR1$variable)
DR1$variable <- gsub("Marriage_Exogamous", "Exogamous Marriage", DR1$variable)
DR1$variable <- gsub("D_Place_", "", DR1$variable)
DR1$variable <- gsub(" definite", "", DR1$variable)
DR1$variable <- gsub("island", "Island", DR1$variable)
DR1$variable <- gsub("_Pattern", "", DR1$variable)
DR1$variable <- gsub("EA033_Extended", "Political Complexity", DR1$variable)

# Aggregate DR test results and write summary to CSV
resDR <- aggregate(DR1[, c("AIC", "Parameters")], by = list(DR1$variable), mean)
names(resDR) <- c("Variable", "Difference", "Overlap")
write.csv(resDR, file = here("outputs", "DR_test_summary.csv"))

#######################################
## SECTION 5b: DR Test % Support     ##
#######################################

# "Parameters" in DR_test: 1 = CIs overlap, 0 = non-overlap
DR1$nonoverlap <- ifelse(DR1$Parameters == 0, 1, 0)
DR1$dummy <- 1

# Match DR_test_summary approach but by variable+location
resDR_loc <- aggregate(DR1[, c("AIC", "Parameters")],
                       by = list(DR1$variable, DR1$location), mean)
names(resDR_loc) <- c("Variable", "Location", "Difference", "Overlap")
resDR_loc$Support_percent <- round(100 * (1 - resDR_loc$Overlap), 2)

# Add counts for transparency
counts_nonoverlap <- xtabs(DR1$nonoverlap ~ DR1$variable + DR1$location)
counts_total      <- xtabs(DR1$dummy      ~ DR1$variable + DR1$location)
counts_df <- data.frame(Variable = rownames(counts_total)[row(counts_total)],
                        Location = colnames(counts_total)[col(counts_total)],
                        N_nonoverlap = as.vector(counts_nonoverlap),
                        N_total = as.vector(counts_total),
                        stringsAsFactors = FALSE)

resDR_loc <- merge(resDR_loc, counts_df, by = c("Variable","Location"), all.x = TRUE)

write.csv(resDR_loc,
          file = here("outputs", "bisse_DR_support_by_location.csv"),
          row.names = FALSE)

#######################################
## SECTION 6: Summarise & Plot BISSE  ##
#######################################

# Define BISSE tests to subset
bisse1 <- c("dull_null", "BiSSE")
bisp <- lf3[lf3$name %in% bisse1, ]

# Add parameter names using the loaded object bisse_param_names.r 
bisp$param_names <- c(xxx, NA) 
bisp <- bisp[bisp$param_names %in% c("turnover0A", "turnover1A", "eps0A", "eps1A"), ]

# Create unique identifiers
bisp$uni  <- paste(bisp$name, bisp$param_names, sep = "_")
bisp$uni2 <- paste(bisp$tree, bisp$variable, bisp$location, sep = ";")

# Rename parameter names for clarity
bisp$param_names[bisp$param_names == "turnover1A"] <- "State 1"
bisp$param_names[bisp$param_names == "turnover0A"] <- "State 0"

# AIC results aggregation and write to CSV
resAIC <- aggregate(bisp$AIC, by = list(bisp$variable, bisp$name, bisp$location), mean)
names(resAIC)[ncol(resAIC)] <- "mean_AIC"
write.csv(resAIC, file = here("outputs", "AIC_test_summary.csv"))

#######################################
## SECTION 7: Parameter Comparison  ##
#######################################

# Compare parameters for BiSSE models between State 1 and State 0
param1 <- bisp[bisp$name == "BiSSE" & bisp$param_names == "State 1", ]
param0 <- bisp[bisp$name == "BiSSE" & bisp$param_names == "State 0", ]
identical(param1$uni2, param0$uni2)

# Calculate mean and standard deviation for each state
mean1 <- aggregate(param1$Parameters, by = list(param1$variable, param1$location), mean)
sd1   <- aggregate(param1$Parameters, by = list(param1$variable, param1$location), sd)
mean0 <- aggregate(param0$Parameters, by = list(param0$variable, param1$location), mean)
sd0   <- aggregate(param0$Parameters, by = list(param0$variable, param1$location), sd)

# Create a results dataframe with ratios of state parameters
resParam <- mean1[, c("Group.1", "Group.2")]
resParam$zero_over_one <- mean0$x / mean1$x
resParam$one_over_zero <- mean1$x / mean0$x
write.csv(resParam, file = here("outputs", "bisse_parameter_differences.csv"))

# Determine where State 1 is greater than State 0
param1$diff <- param1$Parameters - param0$Parameters
param1$diff[param1$diff > 0] <- 1
param1$diff[param1$diff < 0] <- 0
param1$dummy <- 1

# Aggregate counts for the proportion of runs with higher State 1 parameters
x <- xtabs(param1$diff ~ param1$variable + param1$location)
y <- xtabs(param1$dummy ~ param1$variable + param1$location)
z <- x / y
write.csv(z, file = here("outputs", "bisse_parameter_directions.csv"))

#######################################
## SECTION 8: Diversification Analysis ##
#######################################

# Reshape data from long to wide format for diversification calculations
bisp_wide <- dcast(bisp, uni2 ~ uni, value.var = "Parameters")

# Calculate diversification for both BISSE and null models
bisp_wide$diversification0A <- (bisp_wide$BiSSE_turnover0A - (bisp_wide$BiSSE_turnover0A * bisp_wide$BiSSE_eps0A)) / (1 + bisp_wide$BiSSE_eps0A)
bisp_wide$diversification1A <- (bisp_wide$BiSSE_turnover1A - (bisp_wide$BiSSE_turnover1A * bisp_wide$BiSSE_eps1A)) / (1 + bisp_wide$BiSSE_eps1A)
bisp_wide$nulldiversification0A <- (bisp_wide$dull_null_turnover0A - (bisp_wide$dull_null_turnover0A * bisp_wide$dull_null_eps0A)) / (1 + bisp_wide$dull_null_eps0A)
bisp_wide$nulldiversification1A <- (bisp_wide$dull_null_turnover1A - (bisp_wide$dull_null_turnover1A * bisp_wide$dull_null_eps1A)) / (1 + bisp_wide$dull_null_eps1A)

# Select relevant columns and reshape to long format
bisp_wide2 <- bisp_wide[, c("uni2", "diversification0A", "diversification1A", "nulldiversification0A", "nulldiversification1A")]
bisp_long3 <- melt(bisp_wide2, id.vars = "uni2")
bisp_long4 <- cbind(bisp_long3, read.table(text = bisp_long3$uni2, sep = ";", stringsAsFactors = FALSE))
bisp_long5 <- cbind(bisp_long4, read.table(text = bisp_long4$V1, sep = "_", stringsAsFactors = FALSE))
names(bisp_long5) <- c("uni2", "param_names", "Parameters", "uni1", "variable", "area", "type", "tree", "area1", "type2")

# Extract location from uni2 column
bispcols1 <- read.table(text = bisp_long5$uni2, sep = ";")
bispcols2 <- read.table(text = bispcols1$V1, sep = "_")
bisp_long5$location <- bispcols2$V3

# Clean variable names for plotting
bisp_long6 <- bisp_long5[!is.na(bisp_long5$Parameters), ]
bisp_long6$variable <- gsub("Marriage_Endogamous", "Endogamous Marriage", bisp_long6$variable)
bisp_long6$variable <- gsub("Marriage_Exogamous", "Exogamous Marriage", bisp_long6$variable)
bisp_long6$variable <- gsub("D_Place_", "", bisp_long6$variable)
bisp_long6$variable <- gsub(" definite", "", bisp_long6$variable)
bisp_long6$variable <- gsub("island", "Island", bisp_long6$variable)
bisp_long6$variable <- gsub("_Pattern", "", bisp_long6$variable)
bisp_long6$variable <- gsub("EA033_Extended", "Political Complexity", bisp_long6$variable)

# Rename parameter names for diversification
bisp_long6$param_names <- gsub("nulldiversification0A", "Constant Rate", as.character(bisp_long6$param_names))
bisp_long6$param_names <- gsub("nulldiversification1A", "Constant Rate", as.character(bisp_long6$param_names))
bisp_long6$param_names <- gsub("diversification0A", "Trait Absent", as.character(bisp_long6$param_names))
bisp_long6$param_names <- gsub("diversification1A", "Trait Present", as.character(bisp_long6$param_names))
bisp_long6$Parameter <- bisp_long6$param_names

# Create unique identifier and remove duplicates
bisp_long6$uni3 <- paste0(bisp_long6$uni2, bisp_long6$Parameter)
bisp_long6 <- bisp_long6[!duplicated(bisp_long6$uni3), ]

# Quick histogram of Parameters (for inspection)
hist(bisp_long6$Parameters)

#######################################
## SECTION 9: Comparison to NULL Model ##
#######################################

# Extract minimal AIC per tree-variable-location for each model to avoid many-to-many joins
lfB <- lf3[lf3$name == "BiSSE", c("tree","variable","location","AIC")]
lfBn <- lf3[lf3$name == "dull_null", c("tree","variable","location","AIC")]

# One row per key with the minimal AIC (AIC is replicated across params, so min or first is fine)
lfB_min  <- aggregate(AIC ~ tree + variable + location, data = lfB,  FUN = min)
lfBn_min <- aggregate(AIC ~ tree + variable + location, data = lfBn, FUN = min)
names(lfB_min)[names(lfB_min)=="AIC"]   <- "AIC_BiSSE"
names(lfBn_min)[names(lfBn_min)=="AIC"] <- "AIC_null"

# Safe 1-1 merge on all keys
bisse2 <- merge(lfB_min, lfBn_min, by = c("tree","variable","location"), all.x = TRUE, all.y = FALSE)

# Compute delta and flags
bisse2$delta_AIC <- bisse2$AIC_BiSSE - bisse2$AIC_null
bisse2$lessAIC <- as.integer(bisse2$delta_AIC <= -2)
bisse2$dummy <- 1

# Aggregate proportions by variable and location
res2  <- aggregate(bisse2$lessAIC, by = list(bisse2$variable, bisse2$location), sum)
res2a <- aggregate(bisse2$dummy,   by = list(bisse2$variable, bisse2$location), sum)
res2$prop1 <- res2$x / res2a$x
names(res2) <- c("Variable", "Location", "Number", "Proportion")
write.csv(res2, file = here("outputs", "bisse_null_AIC_proportions.csv"))

#######################################
## SECTION 10: LOO Ratios + Support  ##
#######################################

# Merge parameter ratios (by Variable, Location) with DR-test Support for a single deliverable
if (exists("resParam") && exists("resDR_loc")) {
  resParam2 <- resParam
  names(resParam2)[1:2] <- c("Variable","Location")
  # Clean Variable names to match DR tables
  resParam2$Variable <- gsub("Marriage_Exogamous", "Exogamous Marriage", resParam2$Variable)
  resParam2$Variable <- gsub("EA033_Extended", "Political Complexity", resParam2$Variable)
  resParam2$Variable <- gsub("D_Place_", "", resParam2$Variable)
  resParam2$Variable <- gsub(" definite", "", resParam2$Variable)
  resParam2$Variable <- gsub("_Pattern", "", resParam2$Variable)

  loo_numbers <- merge(resParam2,
                       resDR_loc[, c("Variable","Location","Support_percent","N_nonoverlap","N_total")],
                       by = c("Variable","Location"), all.x = TRUE)

  # Rename columns for clarity
  names(loo_numbers)[names(loo_numbers)=="zero_over_one"] <- "ratio_state0_over_state1"
  names(loo_numbers)[names(loo_numbers)=="one_over_zero"] <- "ratio_state1_over_state0"

  write.csv(loo_numbers,
            file = here("outputs", "bisse_LOO_ratio_and_support.csv"),
            row.names = FALSE)
}
