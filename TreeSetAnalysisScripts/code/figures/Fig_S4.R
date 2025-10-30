# -----------------------------------------------------------------------------|

library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(here)

# -----------------------------------------------------------------------------|
# 1) Slopes across models removing specific regions ----
results_folder = "regression_results_regions_fixed"

# List all files matching "regressions25" pattern
files <- list.files(here(results_folder), pattern="regressions25", full.name=TRUE)
files2 <- list.files(here(results_folder), pattern="regressions25", full.name=TRUE)

my_levels <- c("Human Pop. Density","Landscape Friction", "Mean Spoken Area",
               "Distance to City", "Prop. Island Languages")

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

# Filter results to model "include_REGION" only
res1 <- res1[res1$model=="include_REGION",]

# Define model type for convenience
res1$Model <- res1$type

# Create a unique identifier
res1$uni <- paste(res1$tree_type, 
                  res1$type, 
                  res1$tree, 
                  res1$loos_type, 
                  res1$model, 
                  res1$model2, 
                  sep="_")

# Find minimum AIC per unique grouping
res2 <- aggregate(res1$AIC,
                  by=list(res1$tree_type,res1$type,res1$tree,res1$loos_type),
                  FUN=min, na.rm=TRUE)

names(res2)[1:5] <- c("tree_type","type","tree","loos","min")
res2$uni <- paste0(res2$tree_type, res2$type, res2$tree, res2$loos_type)

# Merge min AIC back to main results
res3 <- merge(res1, res2[, c("uni","min")], by="uni", all.x=TRUE)
res3$deltaAIC <- res3$AIC - res3$min

# Determine whether parameter estimates cross zero
res3$p.value <- 1
res3$p.value[res3$X0.025quant<0 & res3$X0.5quant<0 & res3$X0.975quant<0 ] <- 0
res3$p.value[res3$X0.025quant>0 & res3$X0.5quant>0 & res3$X0.975quant>0 ] <- 0

# Make year a factor
res3$Year <- as.factor(res3$Time2)

# Clean variable names
res3$name <- gsub("I\\(","",res3$name,fixed=FALSE)
res3$name <- gsub("\\)","",res3$name,fixed=FALSE)
res3$name <- gsub("\\(","",res3$name,fixed=FALSE)
res3$name <- gsub("1/","",res3$name,fixed=TRUE)
res3$name <- gsub("area_log","Mean Spoken Area",res3$name,fixed=TRUE)
res3$name <- gsub("island","Prop. Island Languages",res3$name,fixed=TRUE)
res3$name <- gsub("distancetocityyear2","Distance to City",res3$name,fixed=TRUE)
res3$name <- gsub("friction_log","Landscape Friction",res3$name,fixed=TRUE)
res3$name <- gsub("popd_log","Human Pop. Density",res3$name,fixed=TRUE)

res3x <- res3[res3$name!="clade_age" & res3$name!="clade_age2" & res3$name!="Intercept" & res3$tree_type=="tree",]

# Factor levels for plotting
res3x$name <- factor(res3x$name,
                     levels=my_levels
)
# -----------------------------------------------------------------------------|
# 2) Slopes across models without kfolding  ----
# -----------------------------------------------------------------------------|
results_folder = "regression_results_fin_log"

# List all files matching "regressions25" pattern
files <- list.files(here(results_folder), pattern="regressions25", full.name=TRUE)
files2 <- list.files(here(results_folder), pattern="regressions25", full.name=TRUE)

my_levels <- c("Human Pop. Density","Landscape Friction", "Mean Spoken Area",
               "Distance to City", "Prop. Island Languages")

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

# Filter results to model "include_REGION" only
res1 <- res1[res1$model=="include_REGION",]

# Define model type for convenience
res1$Model <- res1$type

# Create a unique identifier
res1$uni <- paste(res1$tree_type, 
                  res1$type, 
                  res1$tree, 
                  res1$loos_type, 
                  res1$model, 
                  res1$model2, 
                  sep="_")

# Find minimum AIC per unique grouping
res2 <- aggregate(res1$AIC,
                  by=list(res1$tree_type,res1$type,res1$tree,res1$loos_type),
                  FUN=min, na.rm=TRUE)

names(res2)[1:5] <- c("tree_type","type","tree","loos","min")
res2$uni <- paste0(res2$tree_type, res2$type, res2$tree, res2$loos_type)

# Merge min AIC back to main results
res3 <- merge(res1, res2[, c("uni","min")], by="uni", all.x=TRUE)
res3$deltaAIC <- res3$AIC - res3$min

# Determine whether parameter estimates cross zero
res3$p.value <- 1
res3$p.value[res3$X0.025quant<0 & res3$X0.5quant<0 & res3$X0.975quant<0 ] <- 0
res3$p.value[res3$X0.025quant>0 & res3$X0.5quant>0 & res3$X0.975quant>0 ] <- 0

# Clean variable names in slopes
res1s$name <- gsub("I\\(","",res1s$col,fixed=FALSE)
res1s$name <- gsub("\\)","",res1s$name,fixed=FALSE)
res1s$name <- gsub("\\(","",res1s$name,fixed=FALSE)
res1s$name <- gsub("1/","",res1s$name,fixed=TRUE)
res1s$name <- gsub("area_log","Mean Spoken Area",res1s$name,fixed=TRUE)
res1s$name <- gsub("island","Prop. Island Languages",res1s$name,fixed=TRUE)
res1s$name <- gsub("distancetocityyear2","Distance to City",res1s$name,fixed=TRUE)
res1s$name <- gsub("friction_log","Landscape Friction",res1s$name,fixed=TRUE)
res1s$name <- gsub("popd_log","Human Pop. Density",res1s$name,fixed=TRUE)

# Make year a factor
res3$Year <- as.factor(res3$Time2)

# Subset to tree-based results (excluding clade_age)
res3$name <- gsub("I\\(","",res3$name,fixed=FALSE)
res3$name <- gsub("\\)","",res3$name,fixed=FALSE)
res3$name <- gsub("\\(","",res3$name,fixed=FALSE)
res3$name <- gsub("1/","",res3$name,fixed=TRUE)
res3$name <- gsub("area_log","Mean Spoken Area", res3$name,fixed=TRUE)
res3$name <- gsub("island","Prop. Island Languages", res3$name,fixed=TRUE)
res3$name <- gsub("distancetocityyear2","Distance to City", res3$name,fixed=TRUE)
res3$name <- gsub("friction_log","Landscape Friction", res3$name,fixed=TRUE)
res3$name <- gsub("popd_log","Human Pop. Density", res3$name,fixed=TRUE)

res3all <- res3[res3$name!="clade_age" & res3$name!="clade_age2" & res3$name!="Intercept" & res3$tree_type=="tree",]

# Factor levels for plotting
res3all$name <- factor(res3all$name,
                     levels=my_levels
)

res3all$kfold = "None"
res3all$kfold = as.factor(res3all$kfold)
res3all$Area = "Global"
# -----------------------------------------------------------------------------|
# 3) Figure plotting --------
# -----------------------------------------------------------------------------|

# Modify name factor for second plot
my_levels <- c("Human Pop. Density","Landscape Friction", "Mean Spoken Area",
               "Distance to City", "Prop. Island Languages")
res3x$name <- factor(res3x$name, levels = my_levels)
res3x$kfold <- as.factor(res3x$kfold)

table(res3x$kfold)

res3x$Area = "Global"
res_combined = rbind(res3all, res3x)

# Define consistent palette (matching Ext Data Fig 4), Americas set to yellow-green
cols = c("None" = "black",
         "Eurasia" = '#e7298a',
         "Africa" = '#7570b3',
         "Oceania" = '#d95f02',
         "Americas" = '#66a11e')
# Ensure factor order matches palette names
res_combined$kfold <- factor(res_combined$kfold, levels = names(cols))

dists <- ggplot(
  res_combined[res_combined$Year %in% c(-4250), ],
  aes(x = mean)
) +
  # fill-only
  # outlines matching fill colours
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +
  geom_density(aes(colour = factor(kfold), group = factor(kfold)),
               fill = NA, linewidth = 0.2, alpha = 0.2, show.legend = FALSE) +
  geom_density(aes(fill = factor(kfold), alpha = factor(kfold)), colour = NA) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_cowplot(10) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols, guide = "none") +  # match fills, no extra legend
  scale_alpha_manual(values = c("None" = 0.9, "Africa" = 0.75, "Eurasia" = 0.65, "Oceania" = 0.55, "Americas" = 0.5)) +
  xlab("Slope Parameter Estimate") + ylab("Frequency") +
  labs(fill = "Macroregion excluded", alpha = "Macroregion excluded") +
  theme(legend.position = "top",
        legend.title = element_text(size = 11),
        legend.margin = margin(10,0,0,0))

plot(dists)
ggsave(here("outputs", "Fig_S4", "FigS4_4250_fixed_20Sep2025.png"), plot = dists,
       width = 5.5, height = 6, units = "in", dpi = 600,
       scale = 1.1)

ggsave(here("outputs", "Fig_S4", "FigS4_4250_fixed_20Sep2025.pdf"), plot = dists,
       width = 5.5, height = 6, units = "in", dpi = 600,
       scale = 1.1)