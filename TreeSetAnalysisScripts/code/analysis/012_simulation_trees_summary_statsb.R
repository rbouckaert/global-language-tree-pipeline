# Script: 014_simulation_trees_summary_statsb.R
# Description: Computes summary statistics for simulation trees and outputs per-tree results.
#------------------------------------------------------------------------------
# Input data:
# - simulationtrees/ [Directory: Simulated tree files, pattern "extant"]

# Output data:
# - summary_sim/ [Directory: Per-tree summary statistics, pattern "summary_*_*.csv"]
#------------------------------------------------------------------------------

# Load required libraries
library(moments)
library(apTreeshape)
library(ape)
library(caper)
library(RPANDA)
library(phangorn)
library(phytools)
library(doParallel)
library(picante)
library(data.table)
library(here) # Added here package for path management

# Read in all the trees
tt <- list.files(here("simulationtrees"), pattern="extant", full.names=TRUE)
tt <- tt[sample(1:length(tt), length(tt))]
tt2 <- gsub("lexicon_families", "lexiconfamilies", tt)
tt2 <- read.table(text=tt2, sep="_")
tt2$V6 <- gsub(".tre", "", tt2$V6)

# Loop over trees --------------------------------------------------------------
for (ii in 1:length(tt)){

    # Read in sample tree
    tree <- ape::read.tree(tt[ii])
    
    # If isnull
    if(is.null(tree)) {break}
    
    # Make results data.frame
    res1 <- data.frame(tree=ii, group=tt2[ii,"V4"], type=tt2[ii,"V5"], num=tt2[ii,"V6"])
  
    # Check if done
    if(paste0("summary_", res1$group, "_", res1$num, "_", res1$type, ".csv") %in% 
       list.files(here("summary_sim"), pattern="summary", full.names=FALSE)) {next}
    
    # Tree length
    res1$tree_length <- sum(tree$edge.length, na.rm=TRUE)
  
    # Crown age
    res1$crown_age <- max(picante::node.age(tree)$ages)
  
    # Crown age in years
    res1$crown_age2 <- res1$crown_age * -1000
  
    # Tip
    res1$richness <- length(tree$tip.label)
  
    # Gamma
    res1$gammma <- ape::gammaStat(tree)
    
    # Balance
    res1$balance <- apTreeshape::colless(as.treeshape(tree))
    
    # Tip level DR
    DR <- picante::evol.distinct(tree, type="equal.splits")$w
    
    # Harmonic mean DR
    res1$HmeanDR <- (1 / mean(1 / DR))
    
    # Mean DR
    res1$meanDR <- mean(log(DR))
    
    # SD DR
    res1$sdDR <- sd(log(DR))
    
    # Median DR
    res1$medianDR <- median(DR)
    
    # DR skew
    res1$DRskew <- moments::skewness(DR)
    
    # DR kurtosis
    res1$DRkurtosis <- moments::kurtosis(DR)
    
    e <- simpleError("test error")
    
    # Find rate constant best
    fit.bd <- res <- tryCatch(suppressWarnings(ape::birthdeath(tree)), error=function(e) e)
  
    if(class(res)[1] == "simpleError"){
      res1$birth1 <- NA
      res1$death1 <- NA
    } else {
      res1$birth1 <- bd(fit.bd)[1]
      res1$death1 <- bd(fit.bd)[2]
    }
  
    # Write results
    write.csv(res1, file=here("summary_sim", paste0("summary_", res1$group, "_", res1$num, "_", res1$type, ".csv")))
  
} # End ii trees

