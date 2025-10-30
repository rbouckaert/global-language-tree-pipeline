# List of files and folders read or loaded in this script:
# 1. langa.shp file in the input_data/langa folder
# 2. Files in the outputs folder with names matching "NEW_EDGE_top_100"
# 3. Files in the outputs folder with names matching "ALL_EDGE_top_100"
# 4. languages_and_dialects_geoFINALupdates6.csv file in the input_data folder
library(pacman)

p_load(data.table, rworldmap, rgeos, maptools, cleangeo, 
       raster, sf, rgdal, viridis, ggnewscale, ggplot2, 
       ggthemes, cowplot, tidyverse, here, sp, lm.beta, ggbeeswarm,magrittr)

# -----------------------------------------------------------|
# 1.Data preparation -----------------------------------------
# -----------------------------------------------------------|

# Get Isolate data ------------------------------------------
# 1. Load and summarise ED data -----------------------------------------------
all_ed <- read.csv(here("outputs", "EDGEscores", "ALL_ED_wTREES.csv")) %>%
  mutate(glottocode = if_else(glottocode == "osse1243", "iron1242", glottocode)) %>%
  group_by(glottocode) %>%
  summarise(ED_mean = mean(ED), .groups = "drop")

# 2. Load Glottolog 4.0 families ----------------------------------------------
lg <- read.csv(here("outputs", "Table_S2", "languoid40.csv"))

lg = lg %>%
  transmute(
    glottocode   = if_else(id == "osse1243", "iron1242", id),
    family_id_4  = family_id,
    parent_id)  %>%  
  left_join(lg %>%  #attach family names
              select(id, name) %>% 
              rename(family_id_4 = id, family_name = name),by = "family_id_4")

# 3. Combine ED with Glottolog and flag isolates ------------------------------
lang_ed <- all_ed %>%
  left_join(lg, by = "glottocode") %>%
  add_count(family_id_4, name = "fam_size") %>%
  mutate(
    fam_size         = if_else(family_id_4 == "", 0L, fam_size),
    True_Isolate     = as.integer(family_id_4 == ""),
    Effective_Isolate = as.integer(fam_size == 1)
  )

# 4. Load threat data ----------------------------------------------------------
threat <- read.csv(here("input_data", "processed_threat_data_frame.csv")) %>%
  mutate(glottocode = if_else(glottocode == "osse1243", "iron1242", glottocode))

# 5. Load top EDGE languages----------------------------------------------------
edge_top <- read.csv(here("outputs", "EDGEscores", "NEW_EDGE_top_100.csv")) %>%
  rename(glottocode = Group.1) %>%
  mutate(glottocode = if_else(glottocode == "osse1243", "iron1242", glottocode))

# 6. Merge all sources and compute final flags --------------------------------
table_s2 <- edge_top %>%
  left_join(lang_ed, by = "glottocode") %>%
  left_join(
    threat %>% select(glottocode, status, glottolog_macroarea),
    by = "glottocode"
  ) %>%
  mutate(
    glottolog_macroarea = recode(glottolog_macroarea,
                                 "Australia" = "Oceania",
                                 "Papunesia" = "Oceania"),
    Language_Isolate    = as.integer(True_Isolate == 1 | Effective_Isolate == 1)
  ) %>%
  select(
    name, glottocode, ED, EDGE, EDGErank,
    ED_mean, family_name, family_id_4,fam_size, 
    True_Isolate, Effective_Isolate, Language_Isolate, status, glottolog_macroarea
  )

# World Map data -----------
# Get world map
sPDF2 <- getMap()
sPDF2 <- clgeo_Clean(sPDF2)  # Needed to fix up some non-closed polygons 

sPDF2$NAME <- as.vector(sPDF2$NAME)
sPDF2$NAME[is.na(sPDF2$NAME)] <- "XXX"
sPDF<-spTransform(sPDF2,CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
cont4 <- st_as_sf(sPDF[sPDF$NAME != "Antarctica", ])

# Langauge data
langa <- st_read(here("input_data", "langa", "langa.shp"), quiet = TRUE)

# Edge data
edge1 = table_s2
edge1 %<>% relocate(glottocode, .before = name)

#read in links to iso639
res5 <- fread(file = here("input_data", "languages_and_dialects_geoFINALupdates6.csv"))
res5$`DPLACE subsistence` = res5$`DPLACE subsistence` - 1
res5$glottocode[res5$glottocode=="osse1243"] <- "iron1242"

setkey(res5, glottocode) #- Dave's
res6 <- res5[edge1]

# Check how many isolates
nrow(res6[res6$Language_Isolate == 1, ]) / nrow(res6)

# Check classification differences
table(res6$Language_Isolate)

# Create spatial object
res6_out <- res6[!is.na(res6$EDGErank), ]
res6_out <- res6_out[!is.na(res6_out$latitude), ]
coordinates(res6_out) <-  ~ longitude + latitude
projection(res6_out) <- projection(sPDF2)
res6_out <- spTransform(res6_out, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# Rename and create isolate column
res6_out$Isolate <- "Non-Isolate"
res6_out$Isolate[res6_out$Language_Isolate == 1] <- "Isolate"
res6_out$Isolate <- as.factor(res6_out$Isolate)

# Create rank
res6_out$Rank <- res6_out$EDGErank
res6_out$MeanRank <- log(res6_out$Rank)

# Make sf
res6_out2 <- st_as_sf(res6_out)

# Check
table(res6_out$Isolate)


# -----------------------------------------------------------|
# 2. World map plot -----------------------------------------
# -----------------------------------------------------------|

p<-ggplot()+
  geom_sf(data = cont4,)+
  scale_fill_manual(values = c('#7570b3','#a6761d','#e7298a','#66a61e','#d95f02','#e6ab02'))+
  theme_map()
p<- p+ new_scale_fill() +new_scale_color() 
p<- p+
  geom_sf(data = res6_out2,aes(col=MeanRank,shape=Isolate),size=3)+ 
  scale_colour_viridis(option="A",direction=(-1),alpha=1,begin=0.925,end=0)+
  
  scale_shape_manual(values=c(16,17))+
  guides(
    #reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE),
    #reverse size order (higher diameter on top) 
    size = guide_legend(reverse = TRUE))  
print(p)


# -----------------------------------------------------------|
# 3.Threat category vs Log mean ED data ---------------------
# -----------------------------------------------------------|

all_ed = read.csv(here("outputs","EDGEscores", "ALL_ED_wTREES.csv"))
threat = read.csv(here("input_data", "processed_threat_data_frame.csv"))
all_ed$glottocode[all_ed$glottocode=="osse1243"] <- "iron1242"
threat$glottocode[threat$glottocode=="osse1243"] <- "iron1242"
names(all_ed); dim(all_ed)

#calculate the mean ED for each glottocode
all_ed_mean <- all_ed %>%
  group_by(glottocode) %>%
  summarise(mean_ED = mean(ED, na.rm = TRUE))

# View the resulting dataframe
print(all_ed_mean)

# Join the status information to the all_ed_mean dataframe
all_ed_mean <- all_ed_mean %>%
  left_join(threat %>% dplyr::select(glottocode, status), by = "glottocode")

# View the resulting dataframe
head(all_ed_mean)
all_ed_mean$log_mean_ED = log(all_ed_mean$mean_ED)

eds = all_ed_mean %>% filter(status != "extinct")
eds$status = factor(eds$status, levels = c("safe", "vulnerable", "definitely endangered", "severely endangered", "critically endangered"))

# Convert status to numeric (assuming 1 = Safe, 2 = Vulnerable, ..., 5 = Critically Endangered)
eds$status_num <- as.numeric(factor(eds$status))

# -----------------------------------------------------------|
# 4. Inset plot ---------------------
# -----------------------------------------------------------|

# Compute overall mean
overall_mean <- mean(eds$log_mean_ED, na.rm = TRUE)

violin = ggplot(eds, aes(x = status, y = log_mean_ED))+
  geom_violin(color = NA, alpha = 1, fill = "grey90") +  # Violin plot
  geom_quasirandom(size = 0.5, alpha = 0.2, method = "pseudorandom", fill = "thistle", color = "grey60", shape = 21) +
  stat_summary(fun = "mean", geom = "point", color = "navy", size = 1, shape = 21, fill = "navy")+
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "navy", size = 0.5) +  
  labs(x = "Glottolog threat category", y = "Log mean ED") +  # Axis labels
  theme_minimal() +  # Minimal theme
  theme(
    axis.text.x = element_text(angle = 0, size = 12),  # Rotate x-axis labels
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),       # White plot area background
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )+
  scale_color_brewer(palette = "Reds") +  # Set color palette
  scale_x_discrete(labels = c("1","2","3","4","5"))

 print(violin)

# -----------------------------------------------------------|
# 5. Combined plot ---------------------
# -----------------------------------------------------------|
p = p + 
  guides(
    color = guide_colorbar(reverse = FALSE, # Set colorbar direction to horizontal
                           direction = "horizontal", 
                           barwidth = 5.5, 
                           barheight = 0.6,
                           title.position = "top",
                           title = "Log Mean Rank"),
    size = guide_legend(reverse = FALSE) # Reverse size legend order
  ) +  theme(legend.position = c(0, 0.8),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11)) 
print(p)
# Combine main plot with inset using cowplot
final_plot <- ggdraw() +
  draw_plot(p, 0.04, 0, 1, 1, width = 1) +  # Main plot (full screen)
  draw_plot(violin, x = 0.01, y = 0.07, width = 0.3, height = 0.4)  # Inset position & size

# Print final plot
print(final_plot)

# Save final plot
ggsave(file = here("outputs","Fig_4", "Fig4_15Oct2025.pdf"), 
       plot = final_plot, width = 8, height = 4.5, dpi = 1000)

ggsave(file = here("outputs","Fig_4", "Fig4_15Oct2025.jpg"), 
       plot = final_plot, width = 8, height = 4.5, dpi = 1000)


# MODELS ---------------------------
# LINEAR REGRESSION -----------------------------------------
# 1.Threat category vs Log mean ED data ---------------------
all_ed = read.csv(here("outputs","EDGEscores", "ALL_ED_wTREES.csv"))
threat = read.csv(here("input_data", "processed_threat_data_frame.csv"))
all_ed$glottocode[all_ed$glottocode=="osse1243"] <- "iron1242"
threat$glottocode[threat$glottocode=="osse1243"] <- "iron1242"
names(all_ed); dim(all_ed)

# 2. Join threat status into all_ed using glottocode -----

all_ed <- all_ed %>%
  left_join(threat %>% dplyr::select(glottocode, status), by = "glottocode")

# 3. Remove extinct taxa and prepare the status predictor -----
 
all_ed <- all_ed %>% filter(status != "extinct") # Remove rows where status is "extinct"

# Convert status into an ordered factor
all_ed$status <- factor(all_ed$status, 
                        levels = c("safe", "vulnerable", "definitely endangered", 
                                   "severely endangered", "critically endangered"))

# Create a numeric version of the status variable (safe = 1, vulnerable = 2, etc.)
all_ed <- all_ed %>% mutate(status_num = as.numeric(status))

# 4. Run a regression (ED ~ status_num) for each tree -----
# group by tree, nest the data for each tree, and then run the model on each subset.
tree_models <- all_ed %>%
  group_by(tree) %>%
  nest() %>%
  mutate(
    # Fit the regression model for each tree subset
    model = purrr::map(data, ~ lm(log(ED) ~ status_num, data = .x)),
    # Compute standardized coefficients using lm.beta
    model_beta = purrr:::map(model, ~ lm.beta(.x))
  )


# 5. Extract the desired coefficients and p-values from each model -----

results <- tree_models %>%
  mutate(
    # Unstandardized slope for status_num
    unstd_coef = purrr::map_dbl(model, ~ coef(.x)["status_num"]),
    # Standardized slope for status_num (from lm.beta)
    std_coef   = purrr::map_dbl(model_beta, ~ .x$standardized.coefficients["status_num"]),
    # p-value for the status_num coefficient
    p_value    = purrr::map_dbl(model, ~ summary(.x)$coefficients["status_num", "Pr(>|t|)"])
  ) %>%
  dplyr::select(tree, unstd_coef, std_coef, p_value)


write.csv(results, here("outputs","EDGE", "ThreatStatus_vs_LogED.csv"), row.names = FALSE)
write.csv(results, here("outputs","Fig_4", "ThreatStatus_vs_LogED.csv"), row.names = FALSE)

mean(results$unstd_coef)
mean(results$std_coef)
mean(results$p_value)

# ORDINAL REGRESSION -----------------------------------------
# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(MASS)    # for the polr() function
library(here)


# 1. Read in Data and Preprocess -----------------------------
# Read in the ED scores and threat status data
all_ed <- read.csv(here("EDscoresARTUR", "ALL_ED_ARTUR_wTREES.csv"))
threat <- read.csv(here("input_data", "processed_threat_data_frame.csv"))

# Recode glottocode
all_ed$glottocode[all_ed$glottocode == "osse1243"] <- "iron1242"
threat$glottocode[threat$glottocode == "osse1243"] <- "iron1242"

# Check column names and dimensions
print(names(all_ed))
print(dim(all_ed))  # Expecting 1000 unique tree values

# Join threat status into all_ed (using glottocode)
all_ed <- all_ed %>%
  left_join(threat %>% dplyr::select(glottocode, status), by = "glottocode")

# Remove extinct taxa
all_ed <- all_ed %>% filter(status != "extinct")

# Convert threat status into an ordered factor
all_ed$status <- factor(all_ed$status, 
                        levels = c("safe", "vulnerable", "definitely endangered", 
                                   "severely endangered", "critically endangered"),
                        ordered = TRUE)

# 2. Run Ordinal Regression for Each Tree --------------------

# For each tree,  fit an ordinal regression model with:
#    response: status (ordered factor)
#    predictor: log(ED)

# Note: In an ordinal regression, the coefficient for log(ED) tells how the log-odds
# of being in a higher threat category change with a one-unit increase in log(ED).
# Its exponential is the odds ratio.

tree_models <- all_ed %>%
  group_by(tree) %>%
  nest() %>%
  mutate(
    # Fit the ordinal regression model using polr; Hess = TRUE requests the Hessian (for SEs)
    model = purrr::map(data, ~ polr(status ~ log(ED), data = .x, Hess = TRUE))
  )

# 3. Extract Model Results for the log(ED) Predictor ------------

# For each tree-model  extract:
# - The coefficient for log(ED)
# - Its standard error and t-value (from which  compute a p-value)
# - The odds ratio (exp(coefficient))

results <- tree_models %>%
  mutate(
    # Coefficient for log(ED)
    coef_logED = purrr::map_dbl(model, ~ coef(.x)["log(ED)"]),
    # Calculate odds ratio = exp(coefficient)
    odds_ratio = exp(coef_logED),
    # Extract standard error and t-value for log(ED)
    se_logED = purrr::map_dbl(model, ~ summary(.x)$coefficients["log(ED)", "Std. Error"]),
    t_value = purrr::map_dbl(model, ~ summary(.x)$coefficients["log(ED)", "t value"]),
    # Compute p-value using the normal approximation (2-sided test)
    p_value = 2 * (1 - pnorm(abs(t_value))),
    p_value2 = purrr::map_dbl(model, ~ coeftest(.x)["log(ED)", "Pr(>|t|)"])
  ) %>%
  dplyr::select(tree, coef_logED, odds_ratio, se_logED, t_value, p_value, p_value2)

write.csv(results, here("outputs", "ThreatStatus_vs_logED_ordinal.csv"), row.names = FALSE)


# 4. Summarize and Interpret the Results ---------------------

mean_coef    <- mean(results$coef_logED, na.rm = TRUE)
mean_odds    <- mean(results$odds_ratio, na.rm = TRUE)
mean_p_value <- mean(results$p_value2, na.rm = TRUE)
mean_t_value <- mean(results$t_value, na.rm = TRUE)
cat("Mean coefficient for log(ED):", mean_coef, "\n")
cat("Mean odds ratio:", mean_odds, "\n")
cat("Mean p-value:", mean_p_value, "\n")
cat("Mean t-value:", mean_t_value, "\n")