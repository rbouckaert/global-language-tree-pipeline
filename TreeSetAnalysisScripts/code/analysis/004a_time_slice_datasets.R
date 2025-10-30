# ------------------------------------------------------------------------------#
#                              Time Slice Datasets ----                         #    
# Input data:
# - input_data/languages_and_dialects_geoFINALupdates6.csv [CSV: Language metadata]
# - input_data/processed_threat_data_frame.csv [CSV: Threat data]
# - input_data/islands_transformed.R [RData: Island polygons]
# - input_data/langa/langa.shp [Shapefile: Language polygons]
# - input_data/final_env_data.csv [CSV: Environmental data for polygons]
# - input_data/final_env_datapoints.csv [CSV: Environmental data for points]
# - all_trees/ [Directory: Cut trees for each group]
# - simulationtrees/ [Directory: Simulated trees for each group]

# Output data:
# - datasets_and_trees/ [Directory: time slice datasets and trees]
#   - es4_*_time_*_2.csv [CSV: time-sliced dataset per time slice]
#   - es4_*_time_*_2.tre [Tree: backbone trees for each time slice]

# ----------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------#
#                              Libraries & Setup ------
# ------------------------------------------------------------------------------#

# General libraries
library(pacman)
p_load("moments", "data.table","here")

# Phylogenetic analysis libraries
p_load(apTreeshape,caper,dismo,RPANDA,e1071,phangorn,phytools, 
       picante, mice,raster,paleotree,diversitree,ape,rgdal,
       sp,tcltk,nlme,ade4,terra,rworldmap)

# Source supporting functions
source(here("code","analysis","AUX1_time_slice_functions2.R"))

## Constant error message
e <- simpleError("test error")

## Cut d function
cutD <- function(x, n) {
  cut(x,
    breaks = c(quantile(x, probs = seq(0, 1, by = 1 / n), na.rm = T)),
    include.lowest = TRUE
  )
}

# ------------------------------------------------------------------------------#
#                          Data Preparation ------
# ------------------------------------------------------------------------------#

## Load tip data ----
res5 <- fread(file = here("input_data", "languages_and_dialects_geoFINALupdates6.csv"))
tmp = fread(file = here("input_data", "processed_threat_data_frame.csv"))

res5$`DPLACE subsistence` <- res5$`DPLACE subsistence` - 1

## Get world map----
wrld_simpl2 <- getMap()

## Get spatial data -----
rr <- res5
coordinates(rr) <- ~ longitude + latitude
projection(rr) <- projection(wrld_simpl2)
rr2 <- over(rr, wrld_simpl2)

## Combine spatial data with dataset -----
res5 <- cbind(res5, rr2)
res5[, "REGION" := GEO3major]
res5[, "SUBREGION" := IMAGE24]

## Load island data -----
load(file = here("input_data", "islands_transformed.R"))
bi2$dummy <- 1
rr3 <- over(rr, bi2)
rr3$dummy[is.na(rr3$dummy)] <- 0
res5[, island := rr3$dummy]
setkey(res5, iso_final)

## Load language shapefiles -----
library(sf)
langa <- st_read(here("input_data", "langa", "langa.shp"))
library(data.table)
langa_data <- as.data.table(langa)
langa <- print(langa)
cols_chosen <- c("Shp_Lng", "Shap_Ar")

## Aggregate area
agg1 <- langa_data[order(LANG_IS), lapply(.SD, mean), by = .(LANG_IS),
                  .SDcols = cols_chosen]
setkey(agg1, LANG_IS)


## Merge with ethnolog
res5 <- agg1[res5]

# Read in all points dataframe
# Env data from ethnolog polygons
env_data3 <- fread(here("input_data","fixed","final_env_data.csv"))
env_data3 <- env_data3[!is.na(env_data3$Time3), ]
setkey(env_data3, LANG_IS)

# Match with glottocode
reslu <- res5[, c("LANG_IS", "glottocode"), with = FALSE]
reslu <- reslu[!is.na(LANG_IS), ]
setkey(reslu, LANG_IS)

# Remerge
env_data2 <- env_data3[reslu, nomatch = 0]

## Work out what is in above
res_temp <- res5[!glottocode %in% env_data2$glottocode, ]

env_datapA <- fread(here("input_data", "fixed","final_env_datapoints.csv"))
env_datap <- env_datapA[!is.na(env_datapA$Time3), ]
env_datap <- env_datap[LANG_IS %in% res_temp$glottocode, ]
env_datap[, glottocode := LANG_IS]

### Combine ethnologue derived data with glotto centroid point data
env_data3 <- rbind(env_data2, env_datap)
env_data3[, uni := paste0(glottocode, Time3)]
setkey(env_data3, "uni")

## Get point wise distance to city
env_dataB <- env_datapA[, c("LANG_IS", "Time3", "distancetocityyear")]
names(env_dataB)[3] <- "distancetocityyear2"
env_dataB[, uni := paste0(LANG_IS, Time3)]
env_dataB[, c("LANG_IS", "Time3") := NULL]
setkey(env_dataB, "uni")

## add to env_data
env_data <- env_data3[env_dataB, nomatch = 0]
env_data[, uni := NULL]


# Perform a left join to merge L1_Users from tmp to env_data based on glottocode
# Use match to align L1_Users from tmp to env_data
tmp = fread(file = here("input_data", "processed_threat_data_frame.csv"))
env_data$L1_Users <- tmp$L1_Users[match(env_data$glottocode, tmp$glottocode)]

# ------------------------------------------------------------------------------#
#                              Tree Analysis -----                                   
# ------------------------------------------------------------------------------#

### Read in trees
tt0 <- list.files(here("all_trees"), full.names = TRUE)

## Tree name
tt2a <- list.files(here("all_trees"), full.names = FALSE)

## Get rid of underscores
tt2a <- gsub("lexicon_families", "lexiconfamilies", tt2a)
tt2a <- gsub("X", "", tt2a)
tt2a <- gsub(".tree", "", tt2a)
tt3a <- read.table(text = tt2a, sep = "_")
tt3a$V5 <- 1

## Simualated trees
st0 <- list.files(here("simulationtrees"), pattern = "extant", full.names = TRUE)

## Tree name
st2a <- list.files(here("simulationtrees"), pattern = "extant", full.names = FALSE)
st2a <- gsub(".tre", "", st2a)
st2a <- gsub("lexicon_families", "lexiconfamilies", st2a)
st3a <- read.table(text = st2a, sep = "_")

## Combine together
tt <- c(tt0, st0)
tt2 <- c(tt2a, st2a)
tt3 <- rbind(tt3a, st3a)

### Only need global trees
tt <- tt[tt3$V3 == "Global"]
tt2 <- tt2[tt3$V3 == "Global"]
tt3 <- tt3[tt3$V3 == "Global", ]


## Reorder sample
sam1 <- sample(1:(length(tt)), replace = FALSE)

## Combine together and reorder
tt <- tt[sam1]
tt3 <- tt3[sam1, ]
tt2 <- tt2[sam1]

## Run loop ----
for (ii in 1:length(tt)) {
  time_start = Sys.time()
  
  ## Check if done
  if (paste("es4_", tt2[[ii]], "time", "_3500", "_2.csv", sep = "") %in%
   list.files(here("datasets_and_trees_fixed"), full.names = FALSE)) {
    next
  }

  if (tt3[ii, "V1"] == "simulatione") {
    next
  }

  # For simulated
  tree <- ape::read.tree(tt[ii])

  ## Put real names on simulated tree
  if (tt3[ii, "V1"] == "simulatione") {
    tree <- ladderize(tree)
    rtp <- res5
    rtp$sort <- sample(1:nrow(rtp), nrow(rtp), replace = FALSE)
    rtp$glottocode <- rtp$glottocode[order(rtp$REGION, rtp$sort)]
    tree$tip.label <- rtp$glottocode
  }

  ## Get node ages
  node_ages <- picante::node.age(tree)$ages
  tot_time <- max(node_ages)

  ## Create time slice points
  tots1 <- tot_time - 3.5 ## 3500 years ago
  tots2 <- tot_time - 5 ## 5000 years ago
  tots3 <- tot_time - 4.25 ## 4250 years ago
  tots <- c(tots1, tots2, tots3)

  ### Loop through time points
  for (j in tots) {
    ## Cut up subclades
    tr2 <- treeSlice(tree, j, orientation = "tipwards")

    if (length(tr2) < 10) {
      next
    }
    ## Backbone tree
    tr1 <- treeSlice(tree, j, orientation = "rootwards")
    tr1 <- drop.tip(tr1, tr1$tip.label[tr1$tip.label %in% tree$tip.label])

    ## Data frame of results
    es2 <- data.frame(tree = 1:length(tr2))

    ## Mean diversification rate
    hm_es <- function(x) (1 / mean(picante::evol.distinct(x, type = "equal.splits")$w))
    es2$es_hm <- do.call(c, lapply(tr2, hm_es))

    ## Mean diversification rate
    lm_es <- function(x) (mean(log(1 / picante::evol.distinct(x, type = "equal.splits")$w)))
    es2$es_lm <- do.call(c, lapply(tr2, lm_es))

    ## clade size
    len1 <- function(x) Ntip(x)
    es2$length <- do.call(c, lapply(tr2, len1))

    ## add length
    clade_age <- function(x) max(node.age(x)$ages) * -1000
    es2$clade_age <- do.call(c, lapply(tr2, clade_age))

    ## Skew
    sq_es <- function(x) skewness(evol.distinct(x, type = "equal.splits")$w)
    es2$es_sq <- do.call(c, lapply(tr2, sq_es))
    es2$es_sq[is.na(es2$es_sq)] <- 0

    ## Forager languages
    subs1 <- function(x, res5) {
      res6 <- res5[res5$glottocode %in% x$tip.label, ]
      return(aggregate(res6[, c("DPLACE_subsistence", "Forager_language_definite",
                                "Forager_language_possible", "Shap_Ar", "L1_Users",
                                "island", "latitude", "longitude")],
                                by = list(rep(1, nrow(res6))), mean, na.rm = TRUE))
    }
    subs2 <- unlist(do.call(c, lapply(tr2, subs1, res5 = res5)))

    nn <- 10
    es2$DPLACE <- subs2[seq(from = 2, to = length(subs2), by = nn)]
    es2$Forage1 <- subs2[seq(from = 3, to = length(subs2), by = nn)]
    es2$Forage2 <- subs2[seq(from = 4, to = length(subs2), by = nn)]
    es2$area <- subs2[seq(from = 5, to = length(subs2), by = nn)]
    es2$speakers <- subs2[seq(from = 6, to = length(subs2), by = nn)]
    es2$island <- subs2[seq(from = 7, to = length(subs2), by = nn)]
    es2$latitude <- subs2[seq(from = 9, to = length(subs2), by = nn)]
    es2$longitude <- subs2[seq(from = 10, to = length(subs2), by = nn)]

    ## regions
    subs1x <- function(x, res5) {
      res6 <- res5[res5$glottocode %in% x$tip.label, ]
      return(aggregate(res6[, c("REGION", "SUBREGION")], 
             by = list(rep(1, nrow(res6))), modal, na.rm = TRUE))
    }
    subs2 <- unlist(do.call(c, lapply(tr2, subs1x, res5 = res5)))

    es2$REGION <- subs2[seq(from = 2, to = length(subs2), by = 3)]
    es2$SUBREGION <- subs2[seq(from = 3, to = length(subs2), by = 3)]

    ### Get data per clade from inftemp
    get_dataenv <- function(x, env_data, age, tot_time) {
      maxt <- (tot_time - age) * -1000 

      nearest_time <- env_data$Time3[((env_data$Time3 - maxt)^2) == min((env_data$Time3 - maxt)^2)][1]

      env1b <- env_data[glottocode %in% x$tip.label & Time3 > nearest_time, ]

      if (nrow(env1b) == 0) {
        env3 <- data.frame(env1b[1, -c("glottocode", "LANG_IS")])
        env3[1, ] <- NA
        return(env3)
      }

      env2 <- aggregate(env1b[, -c("glottocode", "LANG_IS")], sum, na.rm = TRUE,
                        by = list(env1b$glottocode))



      return(data.frame(aggregate(subset(env2, select = -c(Group.1)), mean,
             na.rm = TRUE, by = list(rep(1, nrow(env2)))))[, -1])
    }

    ##### END FUNCTION

    ## apply env data function
    env_data_per_tree <- do.call(rbind, lapply(tr2, get_dataenv,
                                              env_data = env_data,
                                              age = j, tot_time = tot_time))
    
    ## scale function
    scale2 <- function(x) {
      if (!is.numeric(x)) {
        return(x)
      }

      if (length(na.omit(x)) == 0) {
        return(rep(0, length(x)))
      }

      if (var(na.omit(x)) == 0) {
        return(x)
      }

      return(scale(x))
    }


    ### Bind back to tree set with group.1
    es2a <- cbind(es2, env_data_per_tree[, -1])

    # scale
    e3a <- apply(es2a[, !names(es2a) %in% c("tree", "REGION", "SUBREGION",
                                            "longitude", "latitude", "length")],
                                            2, scale2)

    ## Non scaled variables
    es3b <- cbind(es2a[, c("tree", "REGION", "SUBREGION", "longitude",
                           "latitude", "length")], e3a)

    ### Inmpute missing values
    temp1 <- suppressWarnings(mice(es3b, m = 1, maxit = 50, method = "pmm",
                                   print = FALSE))
    es3 <- complete(temp1, 1)

    ## Row names to be the same as the backbone tree
    es3$tips <- tr1$tip.label
    row.names(es3) <- tr1$tip.label

    ## change to factor
    es3$region <- as.factor(es3$REGION)
    es3$subregion <- as.factor(es3$SUBREGION)

    ## rename
    es4 <- es3

    ## create data frame
    rownames(es4) <- es4$tips
    es4$temperature2 <- es4$temperature^2
    es4$alt12 <- es4$alt1^2
    es4$dummy <- 1
    es4$clade_age <- (max(es4$clade_age) - es4$clade_age) + 1
    es4$clade_age2 <- es4$clade_age^2

    es4$Time <- j
    es4$Crown <- tot_time
    es4$Time2 <- (tot_time - j) * -1000

    es4$tree_type <- tt3[ii, "V1"]
    es4$tree <- tt3[ii, "V2"]
    es4$group <- tt3[ii, "V3"]

    ## write out files
    write.csv(es4, file = here("datasets_and_trees_fixed", paste("es4_", tt2[[ii]],
                  "time_", -1 * es4$Time2[1], "_2.csv", sep = "")),
                  row.names = FALSE)
    write.tree(tr1, file = here("datasets_and_trees_fixed", paste("es4_", tt2[[ii]],
                   "time_", -1 * es4$Time2[1], "_2.tre", sep = "")))
  } ## j ## end of tots
  time_end = Sys.time()
  time_diff = time_end - time_start
  print(paste("Finished ii:", ii, "in", round(time_diff, 2), "seconds."))
} ## end ii treesBackbone tree
