# Script: 006_regression_summary.R
# Description: Loads regression results from .csv files, processes the results, calculates significance, and summaries
#              Produces summary tables and plots, outputs these tables and plots with modified filenames

# Input data:
# - regression_results_no_alt/ [Directory: Contains regression results]
#   - regressions25_*.csv [CSV: Regression result files]
#   - slopes25_*.csv [CSV: Slope summaries]

# Output data:
# - outputs/ [Directory: Summary tables and plots]
#   - Fig_2/ [Directory: Additional figures]
#   - not_zero_table4_linear.csv [CSV: Table of covariate proportions above/below zero]
# - outputs/Draft_Figures/FigureX_Linear.pdf [PDF: Linear slopes figure]
# - outputs/Draft_Figures/Figure3try4_linear.pdf [PDF: Figure 3, linear slopes]
#------------------------------------------------------------------------------

library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(here)

### Load regression results
# List all files matching "regressions25" pattern
files <- list.files(here("regression_results_fin_log"), pattern="regressions25", full.name=TRUE)
files2 <- list.files(here("regression_results_fin_log"), pattern="regressions25", full.name=TRUE)

res1 <- NULL
for(x in seq_along(files)){
  # Read each regression file
  res1a <- read.csv(files[x], stringsAsFactors = FALSE)
  
  # Extract metadata from filename
  real_data <- read.table(text=files2[x], sep="_", stringsAsFactors = FALSE)[, c(4,6,7)]
  names(real_data) <- c("tree2","year2","type2b")
  real_data$year2 <- gsub("XXX","", real_data$year2)
  
  # Combine
  res1a <- cbind(res1a, real_data)
  
  if(is.null(res1)){
    res1 <- res1a
  } else {
    res1 <- rbind(res1, res1a)
  }
}

# Quick check
plot(res1$tree)

### Filter results to model "include_REGION" only
res1 <- res1[res1$model=="include_REGION",]

### Define model type for convenience
res1$Model <- res1$type

### Create a unique identifier
res1$uni <- paste(res1$tree_type, 
                  res1$type, 
                  res1$tree, 
                  res1$loos_type, 
                  res1$model, 
                  res1$model2, 
                  sep="_")

### Find minimum AIC per unique grouping
res2 <- aggregate(res1$AIC,
                  by=list(res1$tree_type,res1$type,res1$tree,res1$loos_type),
                  FUN=min, na.rm=TRUE)

names(res2)[1:5] <- c("tree_type","type","tree","loos","min")
res2$uni <- paste0(res2$tree_type, res2$type, res2$tree, res2$loos_type)

### Merge min AIC back to main results
res3 <- merge(res1, res2[, c("uni","min")], by="uni", all.x=TRUE)
res3$deltaAIC <- res3$AIC - res3$min

### Determine whether parameter estimates cross zero
res3$p.value <- 1
res3$p.value[res3$X0.025quant<0 & res3$X0.5quant<0 & res3$X0.975quant<0 ] <- 0
res3$p.value[res3$X0.025quant>0 & res3$X0.5quant>0 & res3$X0.975quant>0 ] <- 0

# Clean up variable names for plotting
res3$name <- gsub("I\\(","",res3$name,fixed=FALSE)
res3$name <- gsub("\\)","",res3$name,fixed=FALSE)
res3$name <- gsub("\\(","",res3$name,fixed=FALSE)
res3$name <- gsub("1/","",res3$name,fixed=TRUE)
res3$name <- gsub("area_log","Mean Spoken Area",res3$name,fixed=TRUE)
res3$name <- gsub("island","Prop. Island Languages",res3$name,fixed=TRUE)
res3$name <- gsub("Forage1","Foraging Strategy",res3$name,fixed=TRUE)
res3$name <- gsub("written21","Proportion Written",res3$name,fixed=TRUE)
res3$name <- gsub("distancetocityyear2","Distance to City",res3$name,fixed=TRUE)
res3$name <- gsub("friction_log","Landscape Friction",res3$name,fixed=TRUE)
res3$name <- gsub("statesat","No. of States",res3$name,fixed=TRUE)
res3$name <- gsub("alt_varL","Altitude Variability",res3$name,fixed=TRUE)
res3$name <- gsub("popd_log","Human Pop. Density",res3$name,fixed=TRUE)
res3$name <- gsub("croplandL","Proportion Cropland",res3$name,fixed=TRUE)
res3$name <- gsub("crop_resid","Cropland Residuals",res3$name,fixed=TRUE)
res3$name <- gsub("islandC0","Not island",res3$name,fixed=TRUE)
res3$name <- gsub("islandC1","Island",res3$name,fixed=TRUE)

res3$Year <- as.factor(res3$Time2)

### Subset to tree-based results (excluding clade_age)
res3x <- res3[res3$name!="clade_age" & res3$name!="clade_age2" & res3$name!="Intercept" & res3$tree_type=="tree",]

### Assign geographic area labels
res3x$Area <- "Americas"
res3x$Area[res3x$kfold==3] <- "Eurasia"
res3x$Area[res3x$kfold==1] <- "Africa"
res3x$Area[res3x$kfold==2] <- "Oceania"

### Factor levels for plotting
res3x$name <- factor(res3x$name,
                     levels=c("Human Pop. Density",
                              "Landscape Friction",
                              "Mean Spoken Area",
                              "Distance to City",
                              "Prop. Island Languages")
)

# -----------------------------------------------------------------------------|
# sig_prop5_linear -----------------------------------------------------------
# -----------------------------------------------------------------------------|

### Modify name factor for second plot
my_levels <- c("Distance to City", "Human Pop. Density",
               "Landscape Friction", "Mean Spoken Area", "Prop. Island Languages")

res3x$name2 <- gsub(" Quad.","",res3x$name)
res3x$name2 <- factor(res3x$name2, levels = my_levels)

### Significance checks
res3x$uni2 <- paste(res3x$Year, res3x$name, res3x$tree, sep="_")
res3x$sig <- res3x$sig1

uni2b <- unique(res3x$uni2[res3x$sig==1])
ttt <- read.table(text=uni2b, sep="_")


### Filter out cropland models
res4x <- res3x[res3x$model!="cropland",]
res4x$uni3 <- paste0(res4x$Year, res4x$name2, res4x$tree)
res4x <- res4x[order(res4x$sig, decreasing=TRUE),]
res4x <- res4x[!duplicated(res4x$uni3),]
res4x$dummy <- 1

res4x <- aggregate(res4x[, c("sig","dummy")],
                   by=list(res4x$Year,res4x$name2,res4x$tree),
                   max)
names(res4x)[1:3] <- c("Year","name2","tree")

sigs1 <- aggregate(res4x$sig, by=list(res4x$Year,res4x$name2), sum)
sigs2 <- aggregate(res4x$dummy, by=list(res4x$Year,res4x$name2), sum)
sigs1$pval <- sigs1$x/sigs2$x

### Significance proportion
res3x$uni1 <- paste(res3x$Year, res3x$name, res3x$type, sep="_")
uni1b <- unique(res3x$uni1)

sig_prop <- aggregate(res3x$p.value, by=list(res3x$uni1),
                      function(x) length(x[x==0])/length(x))
sig_prop$x_round = paste0(round(100*sig_prop$x, 2),"%")

filenames_signif = "sig_prop5_linear"

if(!file.exists(here("outputs","fixed", paste0(filenames_signif, ".csv")
))){
  write.csv(sig_prop, file=here("outputs","fixed", paste0(filenames_signif, ".csv")
))
} else {
  print("File already exists!")
}

# -----------------------------------------------------------------------------|
# slopes -----------------------------------------------------------
# -----------------------------------------------------------------------------|


### Calculate mean and sd of slopes
slopesm <- aggregate(res3x$mean, by=list(res3x$uni1), mean, na.rm=TRUE)
slopess <- aggregate(res3x$mean, by=list(res3x$uni1), sd, na.rm=TRUE)
slopes1 <- cbind(slopesm, slopess)

filenames_slopes = "slopes2_linear"

if(!file.exists(here("outputs","fixed", 
                paste0(filenames_slopes, 
                ".csv")))){
    write.csv(slopes1, file=here("outputs","fixed", 
              paste0(filenames_slopes, ".csv")
  ))
} else {
  print("File already exists!")
}

# -----------------------------------------------------------------------------|
# % Not Zero -----------------------------------------------------------
# -----------------------------------------------------------------------------|

### Percentage not zero
tab2 <- NULL
for (yy in seq_along(uni1b)){
  res4y <- res3x[res3x$uni1==uni1b[yy],]
  
  all_obs <- nrow(res4y)
  all_neg <- length(res4y$mean[res4y$X0.025quant<0 & res4y$X0.975quant<0])
  all_pos <- length(res4y$mean[res4y$X0.025quant>0 & res4y$X0.975quant>0])
  
  prop_zero <- all_neg/all_obs
  prop_zero2 <- all_pos/all_obs
  
  tab1 <- data.frame(covariate=uni1b[yy],
                     prop_below_zero=prop_zero,
                     prop_above_zero=prop_zero2,
                     stringsAsFactors = FALSE)
  
  if(is.null(tab2)){tab2 <- tab1} else {tab2 <- rbind(tab2, tab1)}
}

tab2$prop_below_zero2 <- round(tab2$prop_below_zero,3)
tab2$prop_above_zero2 <- round(tab2$prop_above_zero,3)


filename_notzero = "not_zero_table4_linear"

if(!file.exists(here("outputs", "fixed",
                     paste0(filename_notzero = "not_zero_table4_linear", 
                            ".csv")))){
  write.csv(tab2, file=here("outputs", "fixed",
                               paste0(filename_notzero = "not_zero_table4_linear", ".csv")
  ))
} else {
  print("File already exists!")
}
