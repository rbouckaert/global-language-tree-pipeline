# ------------------------------------------------------------------------------#
#            Urban Spatial and Environmental Data Processing              ----
# ------------------------------------------------------------------------------#
# Purpose:
# Processes urban spatial data and environmental variables to extract information
# for language data points or polygons.

# Input data:
# - input_data/urbanspatial-hist-urban-pop-3700bc-ad2000-xlsx.csv [CSV: Historical urban population and location data]
# - input_data/languages_and_dialects_geoFINALupdates6.csv [CSV: Language metadata including coordinates (used for points)]
# - input_data/langa/langa.shp [Shapefile: Language polygons (used for polygons)]
# - input_data/temp_data_10000BC/ [Directory: Raster files for temperature]
# - input_data/gcrop/ [Directory: Raster files for cropland]
# - input_data/2020_walking_only_friction_surface/2020_walking_only_friction_surface.geotiff [Raster: Global friction surface]
# - input_data/baseline/zip/ [Directory: Raster files for population count]
# - input_data/languoid.csv [CSV: Glottolog languoid data]
# - input_data/languages_and_dialects_geoNEW.csv [CSV: Additional language/dialect metadata]

# Output data:
# - outputs/ [Directory: processed environmental outputs]
#   - glotto_languages_cites_states3points.csv [CSV: Main output for point data]
#   - glotto_languages_cites_states3.csv [CSV: Main output for polygon data]
#   - mean_frictionpoints.csv [CSV: Friction data for points]
#   - mean_frictionpoly.csv [CSV: Friction data for polygons]
#   - all_tips_by_year_ethno.csv [CSV: Language data combined with environmental variables]
#   - env_by_year_by_fam_area_all_ethno.csv [CSV: Summarized data by family, area, and overall]
# - input_data/ [Directory: updated environmental data]
#   - final_env_datapoints.csv [CSV: Final environmental data for points]
#   - final_env_data.csv [CSV: Final environmental data for polygons]

#------------------------------------------------------------------------------#
#                            Load Libraries and Data                    ----
# ------------------------------------------------------------------------------#

## Load Libraries ----------------
library(pacman)
p_load("data.table", "raster", "rgdal", "sf", "maptools","exactextractr","here",
  "rnaturalearth", "magrittr","dplyr", "rnaturalearthdata", "raster")
p_load("imputeTS")

## Define Parameters ----------------
# Choose between 'points' or 'polys'
type1 = "polys"
#type1 <- "points"

# ------------------------------------------------------------------------------#
#                             Load Spatial Data                           ----
# ------------------------------------------------------------------------------#

## World polygond using rnaturalearth
ws2 <- ne_countries(scale = "medium", returnclass = "sf")
ws2 <- ws2[ws2$name != "Antarctica", ]

## Set Standard Grid
template <- raster(
  nrow = 360,
  ncol = 864,
  ext = extent(-180, 180, -60, 90),
  crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
)

# ------------------------------------------------------------------------------#
#                           Process Cities Data                           ----
# ------------------------------------------------------------------------------#

### read in SEDAC cities data
cites <- read.csv(
  here(
    "input_data",
    "urbanspatial-hist-urban-pop-3700bc-ad2000-xlsx.csv"
  ),
  stringsAsFactors = FALSE
)

## unique yearly time intervals
cyears <- unique(cites$year)
cyears2 <- seq(from = -4000, to = 2000, by = 500)

## Calculate Distance to Major Cities ----------------
## loop through each year and find distance to major city
for (j in 1:(length(cyears2) - 1)) {
  c2 <- cites[cites$year > cyears2[j] & cites$year <= cyears2[j + 1], ]
  t2 <- aggregate(template, 5)
  t2[cellFromXY(t2, c2[, c("Longitude", "Latitude")])] <- 1
  t3 <- distance(t2)
  names(t3) <- paste("distance_to_city_year_", cyears2[j + 1], sep = "")
  plot(t3, main = cyears2[j + 1])

  if (j == 1) {
    t4 <- t3
  } else {
    t4 <- stack(t4, t3)
  }

  print(j)
}

# ------------------------------------------------------------------------------#
#                          Load Language Data                             ----
# ------------------------------------------------------------------------------#

# read in language data as points or polygons
# points
if (type1 == "points") {
  langa <- read.csv(here("input_data", "languages_and_dialects_geoFINALupdates6.csv"),
                    stringsAsFactors = FALSE)
  langa <- langa[!is.na(langa$latitude) & !is.na(langa$longitude), ]
  coordinates(langa) <- ~ longitude + latitude
}

# polygons
### language shapefiles
if (type1 == "polys") {
  langa <- st_read(here("input_data", "langa", "langa.shp")) |>
    st_make_valid() |>              # repair any self-intersections
    st_cast("MULTIPOLYGON")         # ensure polygon type only
  if (is.na(st_crs(langa))) {
    st_crs(langa) <- 4326                  
  }
}

# ------------------------------------------------------------------------------#
#                          Load Environmental Data                        ----
# ------------------------------------------------------------------------------#

## Read in Temperature Data ----------------
temp1 <- stack(list.files(here("input_data", "temp_data_10000BC"), full.names = TRUE))
projection(temp1) <- projection(raster())
temp2 <- (temp1)

## Read in Cropland Data ----------------


## Read in Friction Surface Data ----------------
fric1 <- stack(c(here("input_data", "2020_walking_only_friction_surface",
                      "2020_walking_only_friction_surface.geotiff")))
projection(fric1) <- projection(raster())


## Read in Population Data ----------------
gpop2 <- raster::stack(list.files(here("input_data", "baseline", "zip"), 
                                  pattern = "popd", full.names = TRUE))
projection(gpop2) <- projection(raster())

## Helper Function
# t4 has 12 layers, so the function must return 12 numbers - one for each layer
trim_mean_layer <- function(values, coverage_fractions) {
  ## 'values' is a matrix: 1 column per layer, many rows (cells) per polygon
  apply(
    X   = values,
    MARGIN = 2,                         # loop over columns / layers
    FUN = function(v) mean((v), trim = 0.10, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------#
#                      Extract Environmental Variables                   ----
# ------------------------------------------------------------------------------#
a = read.csv(here("outputs", "glotto_languages_cites_states3.csv"))

## extract data per language
if (class(langa)[1] == "SpatialPolygonsDataFrame") {
  
  fric2 <- terra::rast(fric1) #NOT LOGGED
  
  ## Extract Mean Friction
  mean_friction <- exact_extract(fric2, langa, "mean") |>
    setNames("mean_friction_walking") |>
    as_tibble()
  
  ## Extract Distance To city
  distance_to_city <- exact_extract(t4, langa, trim_mean_layer)
  distance_to_city <- as.data.frame(t(distance_to_city))   # flip 
  colnames(distance_to_city) <- names(t4) # move layer names to columns
  rownames(distance_to_city) <- NULL
  
  ## Extract Temperature
  temp_at_location <- exact_extract(temp2, langa, trim_mean_layer)
  temp_at_location <- as.data.frame(t(temp_at_location)) # flip 
  temp_at_location <- as.data.frame(temp_at_location)
  colnames(temp_at_location) <- names(temp2) # move layer names to columns
  

  # Extract Population
  # Convert once, then reuse
  gpop2_terra <- terra::rast(gpop2)        # RasterStack  ➜  SpatRaster
  popn_at_location <- exact_extract(gpop2_terra, langa, trim_mean_layer)
  popn_at_location <- as.data.frame(t(popn_at_location)) # flip
  
  ## 1.  attach the friction vector -----------------------------
  # (rename if the column came out as "mean_friction_walking")
  langa$friction <- mean_friction[[1]]    
  
  ## 2.  drop geometry so we have an ordinary data.frame --------
  langa_df <- sf::st_drop_geometry(langa)
  
  ## 3.  bind everything together and write out -----------------
  glot2 <- cbind(
    langa_df,
    distance_to_city,
    temp_at_location,
    popn_at_location
  )
  
  write.csv(
    glot2,
    here("outputs", "fixed", "glotto_languages_cites_states3.csv"),
    row.names = FALSE
  )
  
  write.csv(
    cbind(langa_df, friction = langa$friction),
    here("outputs", "fixed", "mean_frictionpoly.csv"),
    row.names = FALSE
  )
  
} else {
  
  ## For Points ----------------
  
  ## Extract Mean Friction
  mean_friction <- extract((fric1$X2020_walking_only_friction_surface), langa) 
  mean_friction <- as.data.frame(mean_friction)
  names(mean_friction) <- "mean_friction"

  ## Extract Distance To city
  distance_to_city <- extract(t4, langa) 
  distance_to_city <- as.data.frame(distance_to_city)
  names(distance_to_city) <- names(t4)

  ## Extract Temperature
  temp_at_location <- extract(temp2, langa) 
  temp_at_location <- as.data.frame(temp_at_location)
  names(temp_at_location) <- names(temp1)


  ## Extract Population
  popn_at_location <- extract(gpop2, langa) 
  popn_at_location <- as.data.frame(popn_at_location)
  names(popn_at_location) <- names(gpop2)

  ## Put time invariant with langa data
  langa@data$friction <- mean_friction$mean_friction

  ## Write Out Mean Friction Data
  write.csv(cbind(langa@data, mean_friction), here("outputs","fixed", "mean_frictionpoints.csv"))

  ## Combine All Data Together
  glot2 <- cbind(langa@data, distance_to_city, 
                 temp_at_location, popn_at_location)
  write.csv(glot2, here("outputs","fixed", "glotto_languages_cites_states3points.csv"))
}

# ------------------------------------------------------------------------------#
#                          Post-Processing Data                           ----
# ------------------------------------------------------------------------------#

## Read in Processed Data ----------------
# Points
if (type1 == "points") {
  glot3 <- fread(here("outputs","fixed", "glotto_languages_cites_states3points.csv"))
}

# Language shapefiles
if (type1 == "polys") {
  glot3 <- fread(here("outputs","fixed", "glotto_languages_cites_states3.csv"))
}

## Identify Time Variable Columns
# works out where first time variable column is
ttt <- (1:length(names(glot3)))[names(glot3) == "distance_to_city_year_.3500"]

## Reshape Data from Wide to Long
DT1 <- melt(glot3,
  id.vars = names(glot3)[1:(ttt - 1)],
  measure.vars = names(glot3)[ttt:ncol(glot3)]
)

## Process Variable Names
DT1$variable <- gsub("grid_data_", "temperature;", DT1$variable)
DT1$variable <- gsub("distance_to_city_year_", "distancetocityyear;xxx", DT1$variable)
DT1$variable <- gsub("popd_", "popd;", DT1$variable)
DT1$variable <- gsub("cropland", "cropland;", DT1$variable)

## Split Variable and Time
DT1[, c("Variable", "Time") := tstrsplit(variable, ";")]

## Adjust Time Formats
DT1$Time <- gsub(".", "-", DT1$Time, fixed = TRUE)
DT1$Time <- gsub("_years_before_present", ";-1", DT1$Time)
DT1$Time <- gsub("xxx", "1;", DT1$Time)
DT1$Time <- gsub("BP", ";-1", DT1$Time)

## Handle AD/BC Indicators
DT1[, Time2 := Time]
DT1$Time2 <- gsub("AD", ";AD", DT1$Time2)
DT1$Time2 <- gsub("BC", ";BC", DT1$Time2)
DT1$Time <- gsub("AD", ";1", DT1$Time) # now convert
DT1$Time <- gsub("BC", ";-1", DT1$Time)

## Calculate Real Time
DT1[, c("time1", "time2") := tstrsplit(Time, ";")]
DT1[, c("time1x", "time2x") := tstrsplit(Time2, ";")]
DT1[, Time3 := as.numeric(time1) * as.numeric(time2)]
print(length(unique(DT1$Time3)))

### write to file
if (type1 == "points") {
  fwrite(DT1, file = here("outputs","fixed", "interim_env_datapoints_points.csv"))
  DT1 <-fread(file = here("outputs","fixed", "interim_env_datapoints_points.csv"))
} else {
  fwrite(DT1, file = here("outputs","fixed", "interim_env_datapoints.csv"))
  DT1 <- fread(file = here("outputs","fixed", "interim_env_datapoints.csv"))
}

## Rename Language Identifier ----------------
DT1[, LANG_IS := glottocode]
print(length(unique(DT1$Time3)))

## Filter Temperature Data ----------------
## this is because there are two types of temperature data - some BC/AD and some BP
DT2 <- DT1[!(Variable == "temperature" & time2x %in% c("BC", "AD")), ]
print(length(unique(DT2$Time3)))

## Reshape Data from Long to Wide
MyVars <- names(DT2)[!names(DT2) %in% c("friction", "Time3", "LANG_IS", "value", "Variable")]
DT2[, (MyVars) := NULL]
print(length(unique(DT2$Time3)))

DT2[Variable == "popd", "Time3"] <- DT2[Variable == "popd", "Time3"] - max(DT2$Time3[DT2$Variable == "popd"])
print(length(unique(DT2$Time3)))

DT2[Variable == "distancetocityyear", "Time3"] <- DT2[Variable == "distancetocityyear", "Time3"] - max(DT2$Time3[DT2$Variable == "distancetocityyear"])
print(length(unique(DT2$Time3)))

## make separate subtable
DT3 <- DT2[!duplicated(LANG_IS), c("LANG_IS","friction")]
setkey(DT3, LANG_IS) 

## Reshape Data Back to Wide Format
dt3 <- dcast(DT2, Time3 + LANG_IS ~ Variable, value.var = "value", fun.aggregate = function(x) mean(x, na.rm = TRUE))
dt3[Time3 %in% -21000:min(dt3[!is.nan(popd), "Time3"]), "popd"] <- 0

## Choose columns to take forward
cols_chosen <- c("Time3", #"cropland", "friction", 
                 "distancetocityyear", "popd", "temperature")
dt4 <- dt3[order(Time3), lapply(.SD, na_interpolation), by = .(LANG_IS), .SDcols = cols_chosen]

## Set Key for Joining
setkey(dt4, LANG_IS)

# Add in altitude etc
dt4<-dt4[DT3]

## Write Final Environmental Data ----------------
# Points
if (type1 == "points") {
  fwrite(dt4, file = here("input_data","fixed", "final_env_datapoints.csv"))
}

# Language Shapefiles
if (type1 == "polys") {
  fwrite(dt4, file = here("input_data","fixed", "final_env_data.csv"))
  dt4 <- fread(here("input_data","fixed", "final_env_data.csv"))
  sum(!unique(dt4$LANG_IS) %in% unique(lang1$iso639P3code))
}

# ------------------------------------------------------------------------------#
#                           Merge with Language Data                      ----
# ------------------------------------------------------------------------------#

## Read in Additional Language Data ----------------
# add glott0
lang1 <- fread(here("input_data", "languoid.csv"), stringsAsFactors = FALSE)
lang1 <- lang1[!is.na(lang1$iso639P3code) & lang1$iso639P3code != "", ]
setkey(lang1, "iso639P3code")

## Set Key for Environment Data 
setkey(dt4, "LANG_IS")

## Join Language Data 
dt5 <- lang1[dt4]

## Add glott2 data
lang2 <- fread(here("input_data","languages_and_dialects_geoNEW.csv"), stringsAsFactors = FALSE)
setkey(lang2, "glottocode")

## Set Key for Joined Data
setkey(dt5, "id")

## Merge 
dt6 <- lang2[dt5]
dt6[, c("name", "level", "latitude", "longitude")] <- NULL
dt6[, id := glottocode]

## Merge with Ethno Data
setkey(dt6, iso639P3code)

fwrite(dt6, file = here("input_data", "fixed", "final_env_data_glot.csv"))

## Read in Language Shapefiles
langa <- st_read(here("input_data","langa", "langa.shp"), "langa")

## Aggregate Area Data
langa_data <- setDT(langa@data)
cols_chosen <- c("Shp_Lng", "Shap_Ar")
agg1 <- langa_data[order(LANG_IS), lapply(.SD, mean), by = .(LANG_IS),
                   .SDcols = cols_chosen]
setkey(agg1, LANG_IS)

## check
length(unique(dt6$id))

## merge with ethnolog
dt6 <- agg1[dt6]

# check
length(unique(dt6$id))

# ------------------------------------------------------------------------------#
#                         Summarize Data by Groups                        ----
# ------------------------------------------------------------------------------#

## Choose columns for summarising
cols_chosen <- c("Shp_Lng", "Shap_Ar", "cropland", "distancetocityyear", "popd", "temperature")

## Summarise by family
dt_fam <- dt6[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = .(family_id, Time3), .SDcols = cols_chosen]
dt_fam[, type := rep("family", nrow(dt_fam))]
names(dt_fam)[1] <- "variable"

## Summarize by Macroarea
dt_area <- dt6[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = .(macroarea, Time3), .SDcols = cols_chosen]
dt_area[, type := rep("macroarea", nrow(dt_area))]
names(dt_area)[1] <- "variable"

## Summarize Overall
dt_all <- dt6[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = .(Time3), .SDcols = cols_chosen]
dt_all[, type := rep("all", nrow(dt_all))]
dt_all[, variable := rep("all", nrow(dt_all))]

## Combine Summarized Data
dt7 <- rbindlist(list(dt_fam, dt_area, dt_all), use.names = TRUE)

# ------------------------------------------------------------------------------#
#                             Write Final Outputs                          ----
# ------------------------------------------------------------------------------#

fwrite(dt6, file = here("outputs", "all_tips_by_year_ethno.csv"))
fwrite(dt7, file = here("outputs", "env_by_year_by_fam_area_all_ethno.csv"))