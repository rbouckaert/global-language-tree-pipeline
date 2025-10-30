# ------------------------------------------------------------------------------| 
# Load necessary libraries ----
# ------------------------------------------------------------------------------|
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)       # For plot_grid / plot_layout
library(here)
library(psych)
library(dplyr)
library(ggplot2)
library(grDevices)

# ------------------------------------------------------------------------------| 
# Define the CRS to use (equal-area Mollweide) ----
# ------------------------------------------------------------------------------|
my_crs <- st_crs("ESRI:54009")
#aggregation = "mean"
aggregation = "harmonic mean"
#aggregation = "median"

cell_aggregation = "mean"
#cell_aggregation = "median"

# ------------------------------------------------------------------------------|
# Load and process EDGE data and language metadata
# ------------------------------------------------------------------------------|
EDGE2 <- read.csv(here("outputs","EDGEscores", "ALL_EDGE_eq.csv"), stringsAsFactors = FALSE)
# Calculate DR (assuming ED is nonzero)
EDGE2$DR <- 1 / EDGE2$ED

# Load country codes and language metadata; adjust file path as needed
res5 <- read.csv(file = here("input_data", "languages_and_dialects_geoFINALupdates6.csv"),
                 stringsAsFactors = FALSE)
# Fix glottocode if needed
res5$glottocode[res5$glottocode == "osse1243"] <- "iron1242"

# Summarise data per language
if(aggregation == "mean"){
  tt <- aggregate(EDGE2[, c("ED","DR")], by = list(EDGE2$Species), mean, na.rm = TRUE)
} else if(aggregation == "harmonic mean"){
  tt <- aggregate(EDGE2[, c("ED","DR")], by = list(EDGE2$Species), 
                  FUN = function(x) harmonic.mean(x))
} else if(aggregation == "median"){
  tt <- aggregate(EDGE2[, c("ED","DR")], by = list(EDGE2$Species), median, na.rm = TRUE)
}

# Merge with the metadata (using the glottocode/species ID as key)
all_names <- merge(tt, res5[, c("glottocode", "name", "iso_final", "country_ids")],
                   by.x = "Group.1", by.y = "glottocode")
tmp <- all_names %>%
  dplyr::select(Group.1, ED, DR, country_ids) %>%
  rename(lang = Group.1)

# ------------------------------------------------------------------------------|
# Load language spatial data and convert to sf in the chosen CRS
# ------------------------------------------------------------------------------|
langa <- read.csv(here("input_data", "languages_and_dialects_geoFINALupdates6.csv"),
                  stringsAsFactors = FALSE)
langa <- langa[!is.na(langa$latitude) & !is.na(langa$longitude), ]

# (Fix glottocode if needed)
langa$glottocode[langa$glottocode == "osse1243"] <- "iron1242"

# Convert the data frame to an sf object using the latitude and longitude columns.
langa_sf <- st_as_sf(langa, coords = c("longitude", "latitude"), crs = 4326)
# Transform to the equal-area CRS
langa_sf <- st_transform(langa_sf, crs = my_crs)

# Join the EDGE summary (tmp) to the language points by matching glottocode
langa_sf <- left_join(langa_sf, tmp %>% dplyr::select(lang, ED, DR),
                      by = c("glottocode" = "lang"))

 
# Create an equal-area grid (using one CRS throughout) -----

countries <- ne_countries(scale = "medium", returnclass = "sf")
countries <- st_transform(countries, crs = my_crs)
world_bbox <- st_bbox(countries)

# MAKE PLOTS ----------------------------------------------------------------
#for (cell_aggregation in c("mean", "median")) {
for (cell_size in c(200000, 250000, 300000, 400000)){
    
    # Create the grid from the bounding box (converted to an sf polygon)
    grid <- st_make_grid(st_as_sfc(world_bbox), cellsize = cell_size, square = TRUE)
    grid_sf <- st_sf(cell_id = seq_along(grid), geometry = grid, crs = my_crs)
    
    
    # Spatial join: assign each language point to a grid cell ---
    
    points_in_cells <- st_join(langa_sf, grid_sf, left = TRUE)
    
    # Summarize the number of languages and mean diversification rate per cell
    if (cell_aggregation == "mean") {
      agg_fun <- function(x) mean(x, na.rm = TRUE)
    } else if (cell_aggregation == "median") {
      agg_fun <- function(x) median(x, na.rm = TRUE)
    } else {
      stop("cell_aggregation must be either 'mean' or 'median'")
    }
    
    cell_summary <- points_in_cells %>%
      dplyr::group_by(cell_id) %>%
      dplyr::summarize(
        n_lang  = n(),
        mean_dr = agg_fun(DR)
      )
    
    # Join the summary back to the grid
    grid_with_data <- st_join(grid_sf, cell_summary, by = "cell_id")
    grid_with_data <- grid_with_data %>% filter(!is.na(n_lang))
    
    
    # Bivariate color scaling with QUANTILE binning
    # 1) First log-transform and min-max scale each metric:
    # Define a min-max scaling function
    min_max_scale <- function(x) {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    }
    
    grid_with_data <- grid_with_data %>%
      mutate(
        scaled_n_lang  = ifelse(n_lang > 0, (log(n_lang)), NA),
        scaled_mean_dr = ifelse(mean_dr > 0, (log(mean_dr)), NA)
      )
    
    # 2) Use quantiles to create equal-frequency bins for each scaled variable:
    n_bins <- 10
    lang_breaks <- quantile(grid_with_data$scaled_n_lang,
                            probs = seq(0, 1, length.out = n_bins + 1),
                            na.rm = TRUE)
    
    dr_breaks   <- quantile(grid_with_data$scaled_mean_dr,
                            probs = seq(0, 1, length.out = n_bins + 1),
                            na.rm = TRUE)
    
    # Nudge any duplicate breaks by a tiny epsilon
    for (i in seq_along(lang_breaks)[-1]) {
      if (lang_breaks[i] <= lang_breaks[i - 1]) {
        lang_breaks[i] <- lang_breaks[i - 1] + .Machine$double.eps
      }
    }
    
    # Now 'lang_breaks' are guaranteed strictly increasing
    grid_with_data <- grid_with_data %>%
      mutate(
        bin_lang = cut(scaled_n_lang,
                       breaks = lang_breaks,
                       include.lowest = TRUE,
                       labels = FALSE),
        bin_dr   = cut(scaled_mean_dr,
                       breaks = dr_breaks,
                       include.lowest = TRUE,
                       labels = FALSE)
      )
    
    
    # Define a function to interpolate between four corner colors
    interpolate_color <- function(x, y, tl, tr, bl, br) {
      tl_rgb <- grDevices::col2rgb(tl) / 255
      tr_rgb <- grDevices::col2rgb(tr) / 255
      bl_rgb <- grDevices::col2rgb(bl) / 255
      br_rgb <- grDevices::col2rgb(br) / 255
      
      # Horizontal interpolation
      top_rgb    <- (1 - x) * tl_rgb + x * tr_rgb
      bottom_rgb <- (1 - x) * bl_rgb + x * br_rgb
      
      # Vertical interpolation
      final_rgb <- (1 - y) * bottom_rgb + y * top_rgb
      
      final_rgb <- pmin(pmax(final_rgb, 0), 1)
      grDevices::rgb(final_rgb[1], final_rgb[2], final_rgb[3])
    }
    
    # Map the bin indices (1..n_bins) to [0..1] for color interpolation
    grid_with_data <- grid_with_data %>%
      rowwise() %>%
      mutate(
        bivariate_color = ifelse(!is.na(bin_lang) & !is.na(bin_dr),
                                 interpolate_color(
                                   x = (bin_lang - 1) / (n_bins - 1),
                                   y = (bin_dr   - 1) / (n_bins - 1),
                                   tl = "#0096EB",   # top-left (deep blue)
                                   tr = "#6E003E",   # top-right (deep maroon)
                                   bl = "#B0B0B0",   # bottom-left (gray)
                                   br = "#FFD700"    # bottom-right (bright yellow)
                                 ),
                                 NA)
      ) %>%
      ungroup()
    
    
    # Plotting: Map, legend, and correlation scatterplot
    # MAP ---------------
    map_plot <- ggplot() +
      geom_sf(data = grid_with_data, aes(fill = bivariate_color), 
              color = NA, size = 0, alpha = 1) +
      geom_sf(data = countries, fill = NA, color = "gray40", size = 0, alpha = 0.5) +
      scale_fill_identity() +
      labs(
        tag = "A."
      ) +
      coord_sf(ylim = c(-6400000, NA)) +
      theme_minimal() +
      theme(plot.tag = element_text(face = "bold", size = 13),
            plot.tag.position = c(0.055, 0.93),axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    # Number of bins in each dimension
    n_bins <- 10
    
    # Midpoints for each tile along 0..1
    tile_midpoints <- seq(0.5 / n_bins, 1 - 0.5 / n_bins, length.out = n_bins)
    
    # Create a data frame of tile centers
    legend_df <- expand.grid(
      x = tile_midpoints,
      y = tile_midpoints
    )
    
    
    legend_df <- legend_df %>%
      rowwise() %>%
      mutate(
        bivariate_color = interpolate_color(
          x = x,
          y = y,
          tl = "#0096EB",
          tr = "#6E003E",
          bl = "#B0B0B0",
          br = "#FFD700"
        )
      ) %>%
      ungroup()
    
    # LEGEND ---------------
    legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = bivariate_color)) +
      geom_tile(alpha = 1, width = 1 / n_bins, height = 1 / n_bins) +
      scale_fill_identity() +
      scale_x_continuous(
        breaks = seq(0, 1, by = 0.1),
        labels = seq(0, 1, by = 0.1)
      ) +
      scale_y_continuous(
        breaks = seq(0, 1, by = 0.1),
        labels = seq(0, 1, by = 0.1)
      ) +
      labs(
        x = "Language Richness (scaled)",
        y = "Diversification Rate (scaled)", 
        tag = "B."
      ) +
      theme_cowplot(10) +
      theme(plot.margin = margin(r = 5, l = 35, b = 10),
            plot.tag = element_text(face = "bold", size = 13),
            plot.tag.position = c(0, 1.05))
    
    # Scatterplot for correlation
    cell_summary_filtered <- grid_with_data %>% filter(!is.na(n_lang))
    correlation_result <- cor.test(
      log(cell_summary_filtered$n_lang),
      log(cell_summary_filtered$mean_dr),
      method = "spearman"
    )
    
    pvalue_raw <- correlation_result$p.value
    
    if (pvalue_raw > 0.05) {
      pvalue <- "> 0.05"
    } else if (pvalue_raw < 2.2e-16) {
      pvalue <- "< 2.2e-16"
    } else if (pvalue_raw < 0.001) {
      pvalue <- "< 0.001"
    } else if (pvalue_raw < 0.01) {
      pvalue <- "< 0.01"
    } else if (pvalue_raw < 0.05) {
      pvalue <- "< 0.05"
    } else {
      pvalue <- round(pvalue_raw, 3)
    }
    
    # COR PLOT ---------------
    cor_plot <- ggplot(cell_summary_filtered, aes(x = log(n_lang)+1, y = log(mean_dr))) +
      geom_point(alpha = 0.3, color = "steelblue4",position = position_jitter(width = 0.4, height = 0.1)) +
      labs(
        x = "Language Richness (log)",
        y = "Diversification Rate (log)",
        tag = "C."
      ) +
      annotate("text",
               x = range(log(cell_summary_filtered$n_lang))[2]*0.75, 
               y = range(log(cell_summary_filtered$mean_dr))[1] +
                 0.15*sum(abs(range(log(cell_summary_filtered$mean_dr)))),
               label = paste("Spearman's rho =", round(correlation_result$estimate, 2),
                             "\np-value", pvalue),
               hjust = 0, vjust = 0, size = 3, fontface = "bold") +
      theme_minimal(10) +
      theme(plot.margin = margin(r = 50, l = 5, b = 10),
            plot.tag = element_text(face = "bold", size = 13),
            plot.tag.position = c(0, 1.045))
    
    bottom_row <- plot_grid(legend_plot, cor_plot, ncol = 2, rel_widths = c(1, 1.2))
    final_plot <- plot_grid(map_plot, bottom_row, ncol = 1, rel_heights = c(1.9, 1))
    
    
    # Save outputHow 
    ggsave(
      here("outputs", "Fig_S7_equalcells","8April_Quantiles_WithoutMinMax",
           paste0("Map_",aggregation,"_per_language_",cell_aggregation,"_per_cell_",
                  format(cell_size, scientific = FALSE),".jpg")), 
      final_plot, width = 9, height = 7.5, units = "in", dpi = 300 #original width was 10, height 8
    )
  }
#}

