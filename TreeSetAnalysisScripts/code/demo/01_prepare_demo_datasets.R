source(here::here("TreeSetAnalysisScripts", "code", "analysis", "path_utils.R"))

# ------------------------------------------------------------------------------#
#                     Demo Time-Slice Dataset Preparation                        #
# ------------------------------------------------------------------------------#
# Purpose:
#   Create a small, demo-friendly subset of the TreeSet time-slice workflow.
#   The script reads three Global trees from demo/all_trees_demo/ and generates
#   time-sliced dataset/tree pairs for the manuscript cut times (3500, 4250,
#   5000 YBP).
#
# Inputs:
# - demo/all_trees_demo/ [Directory: demo Global trees]
# - input_data/languages_and_dialects_geoFINALupdates6.csv
# - input_data/processed_threat_data_frame.csv
# - input_data/islands_transformed.r
# - input_data/langa/langa.shp
# - input_data/final_env_data.csv
# - input_data/final_env_datapoints.csv
#
# Outputs:
# - demo/datasets_and_trees_demo/ [Directory: generated demo CSV/tree pairs]
#   - es4_*_time_*_2.csv
#   - es4_*_time_*_2.tre
# ------------------------------------------------------------------------------#

library(pacman)
pacman::p_load(
  "ape", "data.table", "here", "mice", "moments", "phytools",
  "picante", "rworldmap", "sf", "sp"
)

ts_source("code", "analysis", "AUX1_time_slice_functions2.R")

demo_here <- function(...) {
  here::here("TreeSetAnalysisScripts", "demo", ...)
}

demo_dir_create <- function(...) {
  dir_path <- demo_here(...)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  dir_path
}

set.seed(20260407)

demo_tree_dir <- demo_here("all_trees_demo")
demo_output_dir <- demo_dir_create("datasets_and_trees_demo")

if (!dir.exists(demo_tree_dir)) {
  stop("Demo tree directory is missing: ", demo_tree_dir)
}

tree_files <- list.files(
  demo_tree_dir,
  pattern = "^tree_[0-9]+_Global_all\\.tree$",
  full.names = TRUE
)

tree_files <- sort(tree_files)

if (length(tree_files) == 0) {
  stop("No demo Global trees were found in ", demo_tree_dir)
}

time_offsets <- c("3500" = 3.5, "4250" = 4.25, "5000" = 5.0)

# ------------------------------------------------------------------------------#
#                          Shared Input Preparation                              #
# ------------------------------------------------------------------------------#

res5 <- data.table::fread(
  file = ts_here("input_data", "languages_and_dialects_geoFINALupdates6.csv")
)
tmp <- data.table::fread(
  file = ts_here("input_data", "processed_threat_data_frame.csv")
)

# The demo uses helper functions that expect a few enriched columns to already
# exist on the main language table.
if (!"DPLACE subsistence" %in% names(res5) && "DPLACE_subsistence" %in% names(res5)) {
  res5$`DPLACE subsistence` <- res5$DPLACE_subsistence
}

res5$L1_Users <- tmp$L1_Users[match(res5$glottocode, tmp$glottocode)]
res5$written <- tmp$written[match(res5$glottocode, tmp$glottocode)]

if ("DPLACE subsistence" %in% names(res5)) {
  res5$`DPLACE subsistence` <- res5$`DPLACE subsistence` - 1
}

wrld_simpl2 <- rworldmap::getMap()

rr <- res5
sp::coordinates(rr) <- ~ longitude + latitude
sp::proj4string(rr) <- sp::proj4string(wrld_simpl2)
rr2 <- sp::over(rr, wrld_simpl2)

res5 <- cbind(res5, rr2)
res5[, REGION := GEO3major]
res5[, SUBREGION := IMAGE24]

load(file = ts_here("input_data", "islands_transformed.r"))
bi2$dummy <- 1
rr3 <- sp::over(rr, bi2)
rr3$dummy[is.na(rr3$dummy)] <- 0
res5[, island := rr3$dummy]
data.table::setkey(res5, iso_final)

langa <- sf::st_read(ts_here("input_data", "langa", "langa.shp"), quiet = TRUE)
langa_data <- data.table::as.data.table(langa)

agg1 <- langa_data[
  order(LANG_IS),
  lapply(.SD, mean),
  by = .(LANG_IS),
  .SDcols = c("Shp_Lng", "Shap_Ar")
]
data.table::setkey(agg1, LANG_IS)

res5 <- agg1[res5]

env_data3 <- data.table::fread(ts_here("input_data", "final_env_data.csv"))
env_data3 <- env_data3[!is.na(env_data3$Time3), ]
data.table::setkey(env_data3, LANG_IS)

reslu <- res5[, c("LANG_IS", "glottocode"), with = FALSE]
reslu <- reslu[!is.na(LANG_IS), ]
data.table::setkey(reslu, LANG_IS)

env_data2 <- env_data3[reslu, nomatch = 0]
res_temp <- res5[!glottocode %in% env_data2$glottocode, ]

env_datap_all <- data.table::fread(ts_here("input_data", "final_env_datapoints.csv"))
env_datap <- env_datap_all[!is.na(env_datap_all$Time3), ]
env_datap <- env_datap[LANG_IS %in% res_temp$glottocode, ]
env_datap[, glottocode := LANG_IS]

env_data2 <- dplyr::select(env_data2, -cropland, -statesat, -alt1, -alt_var)
env_data3 <- rbind(env_data2, env_datap)
env_data3[, uni := paste0(glottocode, Time3)]
data.table::setkey(env_data3, "uni")

env_dataB <- env_datap_all[, c("LANG_IS", "Time3", "distancetocityyear")]
names(env_dataB)[3] <- "distancetocityyear2"
env_dataB[, uni := paste0(LANG_IS, Time3)]
env_dataB[, c("LANG_IS", "Time3") := NULL]
data.table::setkey(env_dataB, "uni")

env_data <- env_data3[env_dataB, nomatch = 0]
env_data[, uni := NULL]
env_data$L1_Users <- tmp$L1_Users[match(env_data$glottocode, tmp$glottocode)]

# ------------------------------------------------------------------------------#
#                               Demo Tree Loop                                   #
# ------------------------------------------------------------------------------#

for (tree_file in tree_files) {
  tree_name <- sub("\\.tree$", "", basename(tree_file))
  tree_parts <- utils::read.table(text = tree_name, sep = "_", stringsAsFactors = FALSE)

  tree <- ape::read.tree(tree_file)
  node_ages <- picante::node.age(tree)$ages
  tot_time <- max(node_ages)

  for (time_label in names(time_offsets)) {
    slice_age <- tot_time - time_offsets[[time_label]]

    csv_name <- paste("es4_", tree_name, "time_", time_label, "_2.csv", sep = "")
    tre_name <- paste("es4_", tree_name, "time_", time_label, "_2.tre", sep = "")
    csv_path <- file.path(demo_output_dir, csv_name)
    tre_path <- file.path(demo_output_dir, tre_name)

    if (file.exists(csv_path) && file.exists(tre_path)) {
      next
    }

    tr2 <- phytools::treeSlice(tree, slice_age, orientation = "tipwards")
    if (length(tr2) < 10) {
      next
    }

    tr1 <- phytools::treeSlice(tree, slice_age, orientation = "rootwards")
    tr1 <- ape::drop.tip(tr1, tr1$tip.label[tr1$tip.label %in% tree$tip.label])

    es2 <- data.frame(tree = seq_along(tr2))
    es2$es_hm <- vapply(tr2, hm_es, numeric(1))
    es2$es_lm <- vapply(tr2, lm_es, numeric(1))
    es2$length <- vapply(tr2, len1, numeric(1))
    es2$clade_age <- vapply(tr2, clade_age, numeric(1))
    es2$es_sq <- vapply(tr2, sq_es, numeric(1))
    es2$es_sq[is.na(es2$es_sq)] <- 0

    subs2 <- unlist(lapply(tr2, subs1, res5 = res5))
    nn <- 10
    es2$DPLACE <- subs2[seq(from = 2, to = length(subs2), by = nn)]
    es2$Forage1 <- subs2[seq(from = 3, to = length(subs2), by = nn)]
    es2$Forage2 <- subs2[seq(from = 4, to = length(subs2), by = nn)]
    es2$area <- subs2[seq(from = 5, to = length(subs2), by = nn)]
    es2$speakers <- subs2[seq(from = 6, to = length(subs2), by = nn)]
    es2$island <- subs2[seq(from = 7, to = length(subs2), by = nn)]
    es2$latitude <- subs2[seq(from = 9, to = length(subs2), by = nn)]
    es2$longitude <- subs2[seq(from = 10, to = length(subs2), by = nn)]

    region_data <- unlist(lapply(tr2, subs1x, res5 = res5))
    es2$REGION <- region_data[seq(from = 2, to = length(region_data), by = 3)]
    es2$SUBREGION <- region_data[seq(from = 3, to = length(region_data), by = 3)]

    env_data_per_tree <- do.call(
      rbind,
      lapply(tr2, get_dataenv, env_data = env_data, age = slice_age, tot_time = tot_time)
    )

    es2a <- cbind(es2, env_data_per_tree[, -1])

    # Follow the historical 004a "2" pattern directly:
    # 1) scale the imputation variables
    # 2) run mice on that scaled matrix
    # 3) reconstruct raw/logged variants afterwards for downstream use
    impute_exclusions <- c("tree", "REGION", "SUBREGION", "longitude", "latitude", "length")

    e3a <- apply(
      es2a[, !names(es2a) %in% impute_exclusions],
      2,
      scale2
    )

    es3b <- cbind(
      es2a[, c("tree", "REGION", "SUBREGION", "longitude", "latitude", "length")],
      e3a
    )

    temp1 <- suppressWarnings(
      mice::mice(es3b, m = 1, maxit = 50, method = "pmm", print = FALSE)
    )
    es3 <- mice::complete(temp1, 1)

    restore_from_scaled <- function(original, scaled_values) {
      observed <- original[!is.na(original)]
      if (length(observed) == 0) {
        return(rep(NA_real_, length(scaled_values)))
      }

      observed_var <- stats::var(observed)
      if (is.na(observed_var) || observed_var == 0) {
        return(rep(observed[1], length(scaled_values)))
      }

      scaled_numeric <- as.numeric(scaled_values)
      observed_mean <- mean(observed)
      observed_sd <- stats::sd(observed)
      scaled_numeric * observed_sd + observed_mean
    }

    es3$area_raw <- restore_from_scaled(es2a$area, es3$area)
    es3$popd_raw <- restore_from_scaled(es2a$popd, es3$popd)
    es3$friction_raw <- restore_from_scaled(es2a$friction, es3$friction)
    es3$distancetocityyear2_raw <- restore_from_scaled(
      es2a$distancetocityyear2,
      es3$distancetocityyear2
    )
    es3$island_raw <- restore_from_scaled(es2a$island, es3$island)

    area_log_unscaled <- log1p(pmax(es3$area_raw, 0))
    popd_log_unscaled <- log1p(pmax(es3$popd_raw, 0))
    friction_log_unscaled <- log1p(pmax(es3$friction_raw, 0))

    es3$area_log <- as.numeric(scale2(area_log_unscaled))
    es3$popd_log <- as.numeric(scale2(popd_log_unscaled))
    es3$friction_log <- as.numeric(scale2(friction_log_unscaled))

    es3$tips <- tr1$tip.label
    row.names(es3) <- tr1$tip.label
    es3$region <- as.factor(es3$REGION)
    es3$subregion <- as.factor(es3$SUBREGION)

    es4 <- es3
    rownames(es4) <- es4$tips
    es4$temperature2 <- es4$temperature^2
    es4$dummy <- 1
    es4$clade_age <- (max(es4$clade_age) - es4$clade_age) + 1
    es4$clade_age2 <- es4$clade_age^2

    es4$Time <- slice_age
    es4$Crown <- tot_time
    es4$Time2 <- as.numeric(time_label) * -1

    es4$tree_type <- "tree"
    es4$tree <- tree_parts[1, 2]
    es4$group <- tree_parts[1, 3]

    utils::write.csv(es4, file = csv_path, row.names = FALSE)
    ape::write.tree(tr1, file = tre_path)
  }
}
