# -----------------------------------------------------------------------------|
# This script:
# - Loads regression results from .csv files
# - Processes the results, calculates significance, and summaries
# - Produces summary tables and plots
# - Outputs these tables and plots with modified filenames (appending "_ARTUR")
#
# Requirements:
# - Directory: here("regression_results_no_alt") containing "regressions25" files
# - Input files: e.g., "regressions25_...csv"
# - The final output files (summary tables, figures) are saved in the "outputs" folder
# -----------------------------------------------------------------------------|

pacman::p_load(ggplot2, cowplot, ggpubr, gridExtra, RColorBrewer, reshape2, here, readr, dplyr, magrittr,
scales, ggh4x)

# 1. LOAD FILES -------------------------------------------------------------------
# List all files matching "regressions25" pattern
results_folder1 = here("regression_results_fin_log")
results_folder2 = here("regression_results_correct_scale")

files <- list.files(results_folder1, pattern="regressions25", full.name=TRUE)
files2 <- list.files(results_folder1, pattern="regressions25", full.name=TRUE)

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

# Subset to tree-based results (excluding clade_age)
res3x <- res3[res3$name!="clade_age" & res3$name!="clade_age2" & res3$name!="Intercept" & res3$tree_type=="tree",]

# Assign geographic area labels
res3x$Area <- "Americas"
res3x$Area[res3x$kfold==3] <- "Eurasia"
res3x$Area[res3x$kfold==1] <- "Africa"
res3x$Area[res3x$kfold==2] <- "Oceania"

# Factor levels for plotting
res3x$name <- factor(res3x$name,
                     levels=my_levels)

# Load slopes data (fitted slopes)
files_slopes <- list.files(results_folder2, pattern="slopes25", full.name=TRUE)
res1s <- NULL
res1s <- files_slopes %>%
  lapply(read_csv, show_col_types = FALSE) %>%
  bind_rows()

# Remove clade_age and simulations
res1s <- res1s[res1s$col!="clade_age" & res1s$tree>0,]

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
res1s$name <- gsub("islandC0","Not island",res1s$name,fixed=TRUE)
res1s$name <- gsub("islandC1","Island",res1s$name,fixed=TRUE)
res1s$uni2 <- paste(res1s$yearbp*-1, res1s$name, res1s$tree, sep="_")
res1s$uni <- paste(res1s$yearbp*-1, res1s$name, sep="_")

# Reorder factors for plotting fitted slopes
res1s$name <- factor(res1s$name, levels = my_levels)

# -----------------------------------------------------------------------------|
# 2) Subplot A: Distribution of slope estimates ----

# Modify name factor for second plot
my_levels <- c("Human Pop. Density","Landscape Friction", "Mean Spoken Area",
               "Distance to City", "Prop. Island Languages")
res3x$name <- factor(res3x$name, levels = my_levels)

dists <- ggplot(
  data = res3x[res3x$Year %in% c(-3500, -4250, -5000), ], # all years
  #data = res3x[res3x$Year %in% c(-4250), ], # 4250 Only
  aes(x = mean, fill = factor(Year))
) +
  geom_density(alpha = 0.7, color = NA) +
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_cowplot(11) +
  xlab("Slope Parameter Estimate") +
  ylab("Frequency") +
  ggtitle("A.") +
  # Move or tweak legend as desired:
  theme(legend.position = "none") 

# Custom palette
dists <- set_palette(dists, "aaas") 
plot(dists)

# -----------------------------------------------------------------------------|
# 3) Subplot B: Fitted slopes for all 3 time-points ----
# -----------------------------------------------------------------------------|
# Create unique identifier
res3x$uni2<-paste(res3x$Year,res3x$name,res3x$tree,sep="_")

# Make sig
res3x$sig=res3x$sig1

# Get significant results
uni2b<-unique(res3x$uni2[res3x$sig==1])

# Keep only the years of interest in res1s
res1s_sub <- subset(res1s, yearbp %in% c(3500, 4250, 5000))

# Make yearbp into a factor so facet_grid displays them in columns
res1s_sub$yearbp_fac <- as.factor(res1s_sub$yearbp)
levels(res1s_sub$yearbp_fac)

# Keep levels for color mapping (5000, 4250, 3500)
res1s_sub$color_yearbp_fac <- factor(res1s_sub$yearbp_fac, levels = c(5000, 4250, 3500))


# CHOOSE WHICH SCALE TO USE ----------------------------------------------------
#scale_type = "log_basic" 
#scale_type = "scaled"
#scale_type = "log_adjusted"
scale_type = "raw"

# Trim x at the 99th percentile *per variable* ------------------
trimmed = res1s_sub
trimmed$y = exp(trimmed$y)

# Define custom clipping levels per group
clip_levels <- c(
  "Human Pop. Density" = 0.75,
  "Landscape Friction" = 0.95,
  "Mean Spoken Area" = 0.97,
  "Distance to City" = 0.95,
  "Prop. Island Languages" = 0.995
)

trimmed <- trimmed %>%
  group_by(name) %>%
  mutate(
    clip_threshold = quantile(x, probs = clip_levels[name[1]], na.rm = TRUE),
    x_trim = pmin(x, clip_threshold)
  ) %>%
  ungroup() %>%
  select(-clip_threshold)  # Remove the intermediate column


# Apply scaling transformation if needed
if(scale_type == "scaled") {
  trimmed <- trimmed %>%
    group_by(name) %>%
    mutate(
      x_trim = scale(x_trim)[,1],  # standardize x (mean=0, sd=1)
    ) %>%
    ungroup()
}

# Split after trimming
bg  <- trimmed %>% filter(!uni2 %in% uni2b)
sig <- trimmed %>% filter( uni2 %in% uni2b)

# Re-calculate ranges on the trimmed data
rng_by_name <- trimmed %>%
  group_by(name) %>%
  summarise(
    xmin = min(x_trim, na.rm = TRUE),
    xmax = max(x_trim, na.rm = TRUE),
    .groups = "drop"
  )

lim_for <- function(nm) {
  rng_by_name %>% filter(name == nm) %>% with(c(xmin, xmax))
}

# 2.  Facetted x-scales (each row gets its own limits) -------------

if(scale_type == "raw"){
  x_axis_title = "Covariate on raw scale"
} else if(scale_type == "scaled") {
  x_axis_title = "Covariate (standardized)"
} else {
  x_axis_title = "Covariate on log scale"
}

if (scale_type == "log_basic") {
         x_scales <- facetted_pos_scales(
           x = list(
             name == "Human Pop. Density"     ~ scale_x_continuous(limits = lim_for("Human Pop. Density"),     trans = "log1p", expand = expansion(mult = 0.02)),
             name == "Landscape Friction"     ~ scale_x_continuous(limits = lim_for("Landscape Friction"),     trans = "log1p", expand = expansion(mult = 0.02)),
             name == "Mean Spoken Area"       ~ scale_x_continuous(limits = lim_for("Mean Spoken Area"),       trans = "log1p", expand = expansion(mult = 0.02)),
             name == "Distance to City"       ~ scale_x_continuous(limits = lim_for("Distance to City"),       trans = "log1p", expand = expansion(mult = 0.02)),
             name == "Prop. Island Languages" ~ scale_x_continuous(limits = lim_for("Prop. Island Languages"), trans = "log1p", expand = expansion(mult = 0.02))
           )
         )
         } else if (scale_type == "raw") {
         x_scales <- facetted_pos_scales(
           x = list(
             name == "Human Pop. Density"     ~ scale_x_continuous(limits = lim_for("Human Pop. Density"),    expand = expansion(mult = 0.02),labels  = scales::label_number(scale_cut = scales::cut_short_scale())),
             name == "Landscape Friction"     ~ scale_x_continuous(limits = lim_for("Landscape Friction"),    expand = expansion(mult = 0.02),labels  = scales::label_number(scale_cut = scales::cut_short_scale())),
             name == "Mean Spoken Area"       ~ scale_x_continuous(limits = lim_for("Mean Spoken Area"),      expand = expansion(mult = 0.02),labels  = scales::label_number(scale_cut = scales::cut_short_scale())),
             name == "Distance to City"       ~ scale_x_continuous(limits = lim_for("Distance to City"),      expand = expansion(mult = 0.02),labels  = scales::label_number(scale_cut = scales::cut_short_scale())),
             name == "Prop. Island Languages" ~ scale_x_continuous(limits = lim_for("Prop. Island Languages"), expand = expansion(mult = 0.02),labels  = scales::label_number(scale_cut = scales::cut_short_scale()))
           )
         )
         } else if (scale_type == "log_adjusted") {  
         x_scales <- facetted_pos_scales(
           x = list(
             name == "Human Pop. Density"     ~ scale_x_continuous(
               limits  = lim_for("Human Pop. Density"),
               trans   = "log1p",
               expand  = expansion(mult = .02),
               guide   = guide_axis(n.dodge = 2, check.overlap = TRUE),
               labels  = scales::label_number(scale_cut = scales::cut_short_scale())
             ),
             name == "Landscape Friction"     ~ scale_x_continuous(
               limits  = lim_for("Landscape Friction"),
               trans   = "log1p",
               expand  = expansion(mult = .02),
               guide   = guide_axis(n.dodge = 2, check.overlap = TRUE),
               labels  = scales::label_number(scale_cut = scales::cut_short_scale())
             ),
             name == "Mean Spoken Area"       ~ scale_x_continuous(
               limits  = lim_for("Mean Spoken Area"),
               trans   = "log1p",
               expand  = expansion(mult = .02),
               guide   = guide_axis(n.dodge = 2, check.overlap = TRUE),
               labels  = scales::label_number(scale_cut = scales::cut_short_scale())
             ),
             name == "Distance to City"       ~ scale_x_continuous(
               limits  = lim_for("Distance to City"),
               trans   = "log1p",
               expand  = expansion(mult = .02),
               guide   = guide_axis(n.dodge = 2, check.overlap = TRUE),
               labels  = scales::label_number(scale_cut = scales::cut_short_scale())
             ),
             name == "Prop. Island Languages" ~ scale_x_continuous(
               limits  = lim_for("Prop. Island Languages"),
               trans   = "log1p",
               expand  = expansion(mult = .02),
               guide   = guide_axis(n.dodge = 2, check.overlap = TRUE),
               labels  = scales::label_number(accuracy = .1)   # small range so plain labels fine
             )
           )
         )
         } else if (scale_type == "scaled") {
         x_scales <- facetted_pos_scales(
           x = list(
             name == "Human Pop. Density"     ~ scale_x_continuous(
               limits  = lim_for("Human Pop. Density"),
               expand  = expansion(mult = .02),
               labels  = scales::label_number(accuracy = 0.1)
             ),
             name == "Landscape Friction"     ~ scale_x_continuous(
               limits  = lim_for("Landscape Friction"),
               expand  = expansion(mult = .02),
               labels  = scales::label_number(accuracy = 0.1)
             ),
             name == "Mean Spoken Area"       ~ scale_x_continuous(
               limits  = lim_for("Mean Spoken Area"),
               expand  = expansion(mult = .02),
               labels  = scales::label_number(accuracy = 0.1)
             ),
             name == "Distance to City"       ~ scale_x_continuous(
               limits  = lim_for("Distance to City"),
               expand  = expansion(mult = .02),
               labels  = scales::label_number(accuracy = 0.1)
             ),
             name == "Prop. Island Languages" ~ scale_x_continuous(
               limits  = lim_for("Prop. Island Languages"),
               expand  = expansion(mult = .02),
               labels  = scales::label_number(accuracy = 0.1)
             )
           )
         )
         }

# 3.  Plot ----------------------------------------------------------

fittedvals <- ggplot() +
  geom_line(
    data      = bg,
    aes(x = x_trim, y = y, group = uni2),
    colour    = "lightgrey", # light-grey backdrop: every fitted curve
    linewidth = 0.1,
    alpha     = 0.3
  ) +
  geom_line(# coloured curves: only significant slopes
    data      = sig,
    aes(x = x_trim, y = y, group = uni2, colour = color_yearbp_fac),
    linewidth = 0.1, alpha = 0.3
  ) +
  facet_grid2(
    rows         = vars(name),
    cols         = vars(yearbp_fac),
    scales       = "free",
    independent  = "x",
    axes         = "all"
  ) +
  x_scales +
  #geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  labs(
    x     = x_axis_title,
    y     = if(scale_type == "scaled") "Fitted response (raw scale)" else "Fitted response (raw scale)",
    title = "B."
  ) +
  theme_cowplot(11) +
  theme(
    legend.position      = "none",
    strip.text.y         = element_blank(),
    strip.text.y.right   = element_blank()
  )

fittedvals <- set_palette(fittedvals, "aaas")
plot(fittedvals)


filename_2 <- "FIXED_Fig_2B_CORRECTEDV2"
ggsave(
  filename = here(
    "outputs","Fig_2",
    paste0(filename_2, "_", scale_type ,"_", gsub("-", "_", Sys.Date()), ".jpg")
  ),
  plot   = fittedvals,
  height = 11,
  width  = 12,
  scale = 0.8
)

# 4. Save combined plot --------------------------------------------------------
library(grid)
combined_plot <- grid.arrange(dists, fittedvals, nullGrob(),
                              ncol = 3, widths = c(1, 1.3, 0.01)) #for all years

filename_2 <- "Fig2_rawY"

ggsave(
  filename = here(
    "outputs","Fig_2",
    paste0(filename_2, "_", scale_type ,"_",gsub("-", "_", Sys.Date()), ".png")
  ),
  plot   = combined_plot,
  height = 11,
  width  = 12,
  scale = 1
)
