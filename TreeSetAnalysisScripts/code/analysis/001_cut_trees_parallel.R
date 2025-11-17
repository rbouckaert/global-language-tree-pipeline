# ------------------------------------------------------------------------------#
#               Language Phylogenies Processing and Tree Statistics         ----
# ------------------------------------------------------------------------------#

# Purpose:
# Takes a set of global language phylogenies and cuts them up into continental 
# and taxonomic groups and summarises tree level statistics on each sub tree

# Input data:
# - input_data/global-language-trees-6636-taxa.trees [Nexus file containing posterior trees]
# - input_data/languages_and_dialects_geoFINALupdates6.csv [CSV file containing language metadata]
# - input_data/top27families.csv [CSV file containing top 27 families]
# - input_data/Lexicongraphic_familes.csv [CSV file containing lexicon families]  

# Output data:
# - progress.log [Log: processing progress]
# - edscores/ [Directory: evolutionary distinctiveness scores]
#   - *.csv [CSV: evolutionary distinctiveness scores for each tree]
# - all_trees/ [Directory: cut trees]
#   - *.tree [Tree: cut trees for each group]
# - summary/ [Directory: summary statistics]
#   - *.csv [CSV: summary statistics for each cut tree]
# - outputs/ [Directory: combined outputs]
#   - All_Trees_Summary.csv [CSV: combined summary of all trees]
#   - in_text_data/ [Directory: in-text data]
#     - DR_rates_across_tips.csv [CSV: summary of DR rates]

# ------------------------------------------------------------------------------#
#                            Load Libraries and Data                    ----
# ------------------------------------------------------------------------------#
library(pacman)
p_load("ape", "picante", "moments", "phytools", "here", "data.table",
       "doParallel", "foreach", "caper", "phyloTop", "RPANDA")

# Create necessary directories if they don't exist
for (dir in c("edscores", "all_trees", "summary", "outputs", "outputs/in_text_data")) {
  dir_path <- here(dir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat("Created directory:", dir_path, "\n")
  }
}

# read in all the trees
tt <- read.nexus(here("input_data", "global-language-trees-6636-taxa.trees"))

# read in family and macroareas
res5 <- fread(file = here("input_data",
                          "languages_and_dialects_geoFINALupdates6.csv"))

res5$`DPLACE subsistence` <- res5$`DPLACE subsistence` - 1


# ------------------------------------------------------------------------------#
#                             Process Data                             ----
# ------------------------------------------------------------------------------#

## collate to Oceania ----------------
res5$glottolog_macroarea[res5$glottolog_macroarea == "Australia" |
                           res5$glottolog_macroarea == "Papunesia"] <- "Oceania"

mac1 <- unique(res5$glottolog_macroarea)

for (y in 1:length(mac1)) {
  rowy <- res5[res5$glottolog_macroarea == mac1[y], ]
  res1 <- data.frame(id = y,
                     group = mac1[y],
                     type = "macroarea",
                     tip = rowy$glottocode)
  if (y == 1) {
    mac2 <- res1
  } else {
    mac2 <- rbind(mac2, res1)
  }
  rm(res1, rowy)
}

## read in family data -------------
lang1 <- read.csv(here("input_data", "top27families.csv"),
                  stringsAsFactors = FALSE)


for (y in 1:nrow(lang1)) {
  rowy <- lang1[y, ]
  tips <- strsplit(rowy$tips, ",")[[1]]
  res1 <- data.frame(id = y,
                     group = rowy$Family,
                     type = "top27families", tip = tips)
  if (y == 1) {
    lang2 <- res1
  } else {
    lang2 <- rbind(lang2, res1)
  }
  rm(res1, rowy)
}


## read in lexicon. family data --------
lex1 <- read.csv(here("input_data",
                      "Lexicongraphic_familes.csv"),
                 stringsAsFactors = FALSE)


for (y in 1:nrow(lex1)) {
  ## subset by y
  lexy <- lex1[y, ]
  
  # if not family use parent id
  if (lexy$Glottocode %in% c("sout3152", "semi1276")) {
    tips <- res5$glottocode[res5$parent_id == "sout3152"]
  } else {
    tips <- res5$glottocode[res5$family_id == lexy$Glottocode]
  }
  
  res1 <- data.frame(id = y,
                     group = lexy$Language.Family,
                     type = "lexicon_families",
                     tip = tips)
  if (y == 1) {
    lex2 <- res1
  } else {
    lex2 <- rbind(lex2, res1)
  }
  rm(res1, lexy)
}


## Format data ----------------
all1 <- rbind(mac2, lang2, lex2)
all1 <- rbind(all1, data.frame(id = 999, group = "Global",
                               type = "all", tip = NA))
all1$uni1 <- paste(all1$group, all1$type, sep = "_")

# remove punctuation issues
all1$uni1 <- gsub("Central_Sudanic", "Central-Sudanic", all1$uni1)
all1$uni1 <- gsub("Nuclear_Trans_New_Guinea",
                  "Nuclear-Trans-New-Guinea", all1$uni1)
all1$uni1 <- gsub("Nuclear_Torricelli", "Nuclear-Torricelli", all1$uni1)
all1$uni1 <- gsub("Lower_Sepik-Ramu", "Lower-Sepik-Ramu", all1$uni1)
all1$uni1 <- gsub("Bantu (Atlantic-Congo)",
                  "Bantu-Atlantic-Congo", fixed = TRUE, all1$uni1)
all1$uni1 <- gsub("North America", "North-America", all1$uni1)
all1$uni1 <- gsub("South America", "South-America", all1$uni1)

# ------------------------------------------------------------------------------#
#                              Write Out Data                             ----
# ------------------------------------------------------------------------------#

# Write out data files if they don't exist yet
if (!file.exists(here("input_data", "groups_table.csv"))) {
  write.csv(all1[!duplicated(all1$uni), c("group", "type")],
            file = here("input_data", "groups_table.csv"))
}

grps <- unique(all1$uni1)

# Write out groups data if file doesn't exist
if (!file.exists(here("input_data", "groups_data1.csv"))) {
  write.csv(grps, file = here("input_data", "groups_data1.csv"))
}

# ------------------------------------------------------------------------------#
#                              Parallel Processing                          ----
# ------------------------------------------------------------------------------#

# DEFINE GROUPS TO CUT 
grps <- unique(all1$uni1)

# 1. Define stats calculation function
calculate_tree_stats <- function(tree, ii, grp1) {
  res1 <- data.frame(
    tree = ii,
    group = grp1$group,
    type = grp1$type
  )
  
  DR <- 1 / picante::evol.distinct(tree, type = "equal.splits")$w
  
  res1$tree_length  <- sum(tree$edge.length, na.rm=TRUE)
  res1$crown_age    <- max(picante::node.age(tree)$ages)
  res1$crown_age2   <- res1$crown_age * -1000
  res1$richness     <- length(tree$tip.label)
  res1$gammma       <- ape::gammaStat(tree)
  res1$balance      <- colless.phylo(tree, normalise = FALSE)
  
  # DR statistics
  res1$HmeanDR      <- 1 / mean(1 / DR)
  res1$meanDR       <- mean(DR)
  res1$minDR        <- min(DR)
  res1$maxDR        <- max(DR)
  res1$mean_log_DR  <- mean(log(DR))
  res1$sdDR         <- sd(DR)
  res1$sd_log_DR    <- sd(log(DR))
  res1$medianDR     <- median(DR)
  res1$DRskew       <- moments::skewness(DR)
  res1$DRkurtosis   <- moments::kurtosis(DR)
  
  # Birth-death
  fit.bd <- tryCatch(
    suppressWarnings(ape::birthdeath(tree)),
    error = function(e) e
  )
  if (inherits(fit.bd, "simpleError")) {
    res1$birth1 <- NA
    res1$death1 <- NA
  } else {
    res1$birth1 <- bd(fit.bd)[1]
    res1$death1 <- bd(fit.bd)[2]
  }
  
  return(res1)
}

# 2. Define tree processing function
process_single_tree <- function(ii, tt, grps, all1) {
  
  tree2 <- as(tt[[ii]], "phylo") 
  
  # fix for the Ossetic label
  tree2$tip.label[tree2$tip.label == "osse1243"] <- "iron1242"
  
  ED <- picante::evol.distinct(tree2, type = "equal.splits")
  names(ED) <- c("glottocode", "ED")
  
  ## write out ED
  write.csv(ED, file = here("edscores", paste0("edscores_", ii, ".csv")))
  
  results <- list()
  
  for (jj in seq_along(grps)) {
    # to process each group in 'grps'
    # rather than "Global_all" every time:
    grp1 <- all1[all1$uni1 == grps[jj], ]
    if (nrow(grp1) == 0) next  # skip if no matching rows
    
    # If group is 'all', keep full tree
    if (unique(grp1$type) == "all") {
      tree <- tree2
    } else {
      drop_tips <- setdiff(tree2$tip.label, grp1$tip)
      tree <- ape::drop.tip(tree2, drop_tips)
    }
    
    # Save tree if not present
    ## Check if the file already exists
    file_path <- here("all_trees", paste0("tree_", ii, "_", unique(grp1$uni1), ".tree"))
    
    if (!file_path %in% list.files(here("all_trees"), pattern = "tree", full.names = TRUE)) {
      ape::write.tree(tree, file = file_path)
    }

    summary_path <- here("summary", 
                         paste0("summary_", ii, "_", unique(grp1$uni1), ".csv"))
    
    if (!file.exists(summary_path)) {
      res1 <- calculate_tree_stats(tree, ii, grp1[1, ])
      write.csv(res1, file = summary_path, row.names = FALSE)
      results[[length(results) + 1]] <- res1
    }
  }
  # Progress message -----
  cat(sprintf("Finished iteration %d / 1000\n", ii),
      file = "progress.log", append = TRUE)
  # ------
  return(results)
}

# 3. Setup parallel processing
num_cores <- detectCores() - 2
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 4. Main parallel execution
results <- foreach(
  ii = 1:1000,
  .packages = c("ape", "picante", "moments", "here", "phytools", "phyloTop"),
  .combine = 'c'
) %dopar% {
  process_single_tree(ii, tt, grps, all1)
}

stopCluster(cl)

# 5. Combine results into single data.frame
df_all <- rbindlist(results)

# 6. Save results
write.csv(df_all, file = here("outputs", "All_Trees_Summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------#
#                              Post-processing                              ----
# ------------------------------------------------------------------------------#

tmp = df_all %>% filter(type == "all")

# Function to calculate 95% HDPI using 'coda'
calculate_hdpi_coda <- function(data_vector) {
  # Convert the data vector to an MCMC object
  mcmc_obj <- as.mcmc(data_vector)
  
  # Calculate the HDPI interval
  hdpi <- HPDinterval(mcmc_obj, prob = 0.95)
  
  return(hdpi)
}

# Calculate 95% HDPI for 'medianDR'
hdpi_medianDR <- calculate_hdpi_coda(tmp$medianDR)

# Calculate 95% HDPI for 'minDR'
hdpi_minDR <- calculate_hdpi_coda(tmp$minDR)

# Calculate 95% HDPI for 'maxDR'
hdpi_maxDR <- calculate_hdpi_coda(tmp$maxDR)

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
  Median = c(median(tmp$medianDR, na.rm = TRUE),
             median(tmp$minDR, na.rm = TRUE),
             median(tmp$maxDR, na.rm = TRUE)),
  HDPI_Lower = c(hdpi_medianDR[1],
                 hdpi_minDR[1],
                 hdpi_maxDR[1]),
  HDPI_Upper = c(hdpi_medianDR[2],
                 hdpi_minDR[2],
                 hdpi_maxDR[2])
)

print(summary_table_coda)
summary_table_coda2 = summary_table_coda
summary_table_coda2$Median = round(summary_table_coda2$Median, 2)
summary_table_coda2$HDPI_Lower = round(summary_table_coda2$HDPI_Lower, 2)
summary_table_coda2$HDPI_Upper = round(summary_table_coda2$HDPI_Upper, 2)

write.csv(summary_table_coda, file = here("outputs","in_text_data","DR_rates_across_tips.csv"), row.names = FALSE)
