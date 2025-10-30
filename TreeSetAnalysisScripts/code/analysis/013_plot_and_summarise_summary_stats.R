# Script: 013_plot_and_summarise_summary_stats.R
# Description: Reads summary statistics files and generates summary plots and tables.
#------------------------------------------------------------------------------
# Input data:
# - outputs/ [Directory: Contains summary statistics files, e.g., All_Trees_Summary.csv, summary_statsALL_*, summary_statsSIM_*]
#
# Output data:
# - outputs/ [Directory: Contains summarised results]
#   - summary_statsMEAN_<date>.csv [CSV: Mean summary statistics, date-stamped]
#   - simulated_summary_statsALL_HPDI_<date>.csv [CSV: Simulated summary statistics with HPDI, date-stamped]

#------------------------------------------------------------------------------


# Read in summary files and create plots

library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(ggridges)
library(data.table)
library(forcats)
library(here) 
library(coda)
library(bazar)

# Helper function --------------------
HPDI <- function( samples , prob = 0.95) {
  coerce.list <- c( "numeric" , "matrix" , "data.frame" , "integer" , "array" )
  if ( inherits(samples, coerce.list) ) {
    # Single chain for single variable
    samples <- coda::as.mcmc( samples )
  }
  x <- sapply( prob , function(p) coda::HPDinterval( samples , prob=p ) )
  # Now order inside-out in pairs
  n <- length(prob)
  result <- rep(0,n*2)
  for ( i in 1:n ) {
    low_idx <- n+1-i
    up_idx <- n+i
    # Lower
    result[low_idx] <- x[1,i]
    # Upper
    result[up_idx] <- x[2,i]
    # Add names
    names(result)[low_idx] <- concat("|",prob[i])
    names(result)[up_idx] <- concat(prob[i],"|")
  }
  return(result)
}
# End Helper function

# read in files
# lf1 <- list.files(here("outputs"), pattern="summary_statsALL_2", full.names = TRUE)
# lf2 <- list.files(here("outputs"), pattern="summary_statsALL_2", full.names = FALSE)
# lf2 <- gsub(".csv", "", lf2)
# lf2 <- gsub("ALL_", ";", lf2)
# lf3 <- read.table(text=lf2, sep=";")
# lf3$filen <- lf1
# lf3$V2 <- as.Date(lf3$V2)
# 
# # read in latest file
# fil3 <- fread(lf3[lf3$V2==(sort(lf3$V2, decreasing = TRUE)[1]), "filen"]) # choose the one with the latest date
fil3 = fread(here("outputs","All_Trees_Summary.csv"))
#summarise by group
#mean

tt <- aggregate(fil3[,which(names(fil3) %in% "tree_length"):(ncol(fil3)-2)],
                by=list(fil3$type, fil3$group), mean, na.rm=TRUE)
#sd
tt2 <- aggregate(fil3[,which(names(fil3) %in% "tree_length"):(ncol(fil3)-2)],
                 by=list(fil3$type, fil3$group), sd, na.rm=TRUE)
#median
tt3 <- aggregate(fil3[,which(names(fil3) %in% "tree_length"):(ncol(fil3)-2)],
                 by=list(fil3$type, fil3$group), median, na.rm=TRUE)

nd2 <- NULL
#reorder data frames
for (f in 3:ncol(tt)){
  
  nd <- data.frame(type=tt$Group.1, language_group=tt$Group.2, variable=names(tt)[f], mean=tt[,f], sd=tt2[,f], upper95CI=(tt[,f])+(tt2[,f])/(1000^0.5), lower95CI=(tt[,f])-(tt2[,f])/(1000^0.5), median=tt3[,f])
  
  if(is.null(nd2)){nd2 <- nd} else{nd2 <- rbind(nd2, nd)}
  
}

#add on group
nd3 <- nd2 #merge(nd2, all1, by.x="language_group", by.y="group")

#rename
names(nd3)[1:2] <- c("type", "group")

#sort
nd3 <- nd3[order(nd3$type, nd3$group),]

#swap names
nd3$variable <- gsub("gammma", "Gamma Statistic", nd3$variable)
nd3$variable <- gsub("crown_age2", "Crown Age (Years)", nd3$variable)
nd3$variable <- gsub("HmeanDR", "Harmonic Mean Tip\nDiversification Rate", nd3$variable)
nd3$variable <- gsub("richness", "Number of Languages in Group", nd3$variable)
nd3$variable <- gsub("tree_length", "Phylogenetic Diversity", nd3$variable)
nd3$variable <- gsub("crown_age", "Crown Age (1000 Years Ago)", nd3$variable)
nd3$variable <- gsub("meanDR", "Mean Tip Diversification Rate", nd3$variable)
nd3$variable <- gsub("balance", "Phylogenetic Balance", nd3$variable)
nd3$variable <- gsub("sdDR", "Standard Deviation in\nDiversification Rate", nd3$variable)
nd3$variable <- gsub("medianDR", "Median Tip\nDiversification Rate", nd3$variable)
nd3$variable <- gsub("DRskew", "Skew in\nDiversification Rate ", nd3$variable)
nd3$variable <- gsub("DRkurtosis", "Kurtosis in\nDiversification Rate", nd3$variable)

#write summarised results
write.csv(nd3, file = here("outputs", paste0("summary_statsMEAN_", Sys.Date(), ".csv")))
# nd3 <- fread("./summary_statsMEAN_2023-07-25.csv")

#summarise by group
#mean
# column from which to start calculating
  col_start = which(names(fil3)%in% "tree_length")
tt1 <- aggregate(fil3[,col_start:(ncol(fil3)-2)], by=list(fil3$type, fil3$group), function (x) HPDI(x, 0.95))
#90
tt2 <- aggregate(fil3[,col_start:(ncol(fil3)-2)], by=list(fil3$type, fil3$group), function (x) HPDI(x, 0.90))
#99
tt3 <- aggregate(fil3[,col_start:(ncol(fil3)-2)], by=list(fil3$type, fil3$group), function (x) HPDI(x, 0.99))

#median
tt4 <- aggregate(fil3[,col_start:(ncol(fil3)-2)], by=list(fil3$type, fil3$group), median, na.rm=TRUE)
#median
tt5 <- aggregate(fil3[,col_start:(ncol(fil3)-2)], by=list(fil3$type, fil3$group), min, na.rm=TRUE)
#median
tt6 <- aggregate(fil3[,col_start:(ncol(fil3)-2)], by=list(fil3$type, fil3$group), max, na.rm=TRUE)

nd2 <- NULL
#reorder data frames
for (f in 3:ncol(tt4)){
  
  nd <- data.frame(type=tt1$Group.1, language_group=tt1$Group.2, variable=names(tt1)[(f)], range95=tt1[(f)], range90=tt2[(f)], range99=tt3[(f)], median=tt4[,f], min=tt5[,f], max=tt6[,f])
  names(nd)[4:6] <- c("range95", "range90", "range99")
  
  if(is.null(nd2)){nd2 <- nd} else{nd2 <- rbind(nd2, nd)}
  
}

#add on group
nd3 <- nd2 #merge(nd2, all1, by.x="language_group", by.y="group")

#write results
write.csv(nd3, file = here("outputs", paste0("simulated_summary_statsALL_HPDI_", Sys.Date(), ".csv")))
# nd3 <- fread("./simulated_summary_statsALL_HPDI_2023-07-25.csv")

# read in simulated tress and create plots

#rename
fil3a <- fil3
rm(fil3)

#reorder
fil3a <- fil3a[order(sort(fil3a$meanDR)),]

#fil3$name <- fil3$group
fil3a$type = "real"

#create a single summary file countries for simulated data
lf1 <- list.files(here("outputs"), pattern="summary_statsSIM_2", full.names = TRUE)
lf2 <- list.files(here("outputs"), pattern="summary_statsSIM_2", full.names = FALSE)
lf2 <- gsub(".csv", "", lf2)
lf2 <- gsub("SIM_", ";", lf2)
lf3 <- read.table(text=lf2, sep=";")
lf3$filen <- lf1
lf3$V2 <- as.Date(lf3$V2)

# read in latest file
fil3 <- fread(lf3[lf3$V2==(sort(lf3$V2, decreasing = TRUE)[1]), "filen"]) # choose the one with the latest date

write.csv(fil3, file=paste(here("outputs", "summary_statsALL_", sample(1:100000, 1), ".csv")))

# End of script