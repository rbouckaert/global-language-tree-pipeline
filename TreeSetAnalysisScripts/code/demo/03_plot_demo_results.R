source(here::here("TreeSetAnalysisScripts", "code", "analysis", "path_utils.R"))

# ------------------------------------------------------------------------------#
#                         Demo Results Visualisation                             #
# ------------------------------------------------------------------------------#
# Purpose:
#   Create a compact, demo-friendly summary figure for the TreeSet demo.
#   This is intentionally lightweight: it visualises the small demo regression
#   outputs without attempting to recreate the manuscript figure pipeline.
#
# Caveats:
# - This visualisation is based on 3 demo trees rather than the full posterior.
# - The very high in-sample fit metrics from the demo regressions are not shown
#   as headline results, because they are not meaningful as out-of-sample
#   performance measures.
#
# Inputs:
# - demo/regression_results_demo/*.csv
# - demo/regression_slopes_demo/*.csv
# - demo/demo_summary/demo_regression_run_summary.csv
#
# Outputs:
# - demo/demo_figures/demo_regression_summary.png
# - demo/demo_figures/demo_regression_summary.pdf
# - demo/demo_summary/demo_coefficients_summary.csv
# - demo/demo_summary/demo_session_info.txt
# ------------------------------------------------------------------------------#

pacman::p_load(
  "cowplot", "dplyr", "ggplot2", "here", "readr", "scales", "stringr"
)

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

results_dir <- demo_here("regression_results_demo")
slopes_dir <- demo_here("regression_slopes_demo")
summary_file <- demo_here("demo_summary", "demo_regression_run_summary.csv")
figure_dir <- demo_dir_create("demo_figures")
summary_dir <- demo_dir_create("demo_summary")

result_files <- list.files(results_dir, pattern = "regressions25\\.csv$", full.names = TRUE)
slope_files <- list.files(slopes_dir, pattern = "slopes25\\.csv$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No demo regression result files found in ", results_dir)
}

if (length(slope_files) == 0) {
  stop("No demo slope files found in ", slopes_dir)
}

if (!file.exists(summary_file)) {
  stop("Run summary file is missing: ", summary_file)
}

run_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
if (!all(run_summary$status == "success")) {
  stop("Not all demo regression runs completed successfully; inspect ", summary_file)
}

pretty_names <- c(
  "popd_log" = "Human Pop. Density",
  "friction_log" = "Landscape Friction",
  "area_log" = "Mean Spoken Area",
  "distancetocityyear2" = "Distance to City",
  "island" = "Prop. Island Languages"
)

covariate_levels <- c(
  "Human Pop. Density",
  "Landscape Friction",
  "Mean Spoken Area",
  "Distance to City",
  "Prop. Island Languages"
)

time_levels <- c("3500", "4250", "5000")

coefficients <- bind_rows(lapply(result_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$file <- basename(f)
  d
}))

coefficients <- coefficients %>%
  filter(name %in% names(pretty_names), model == "include_REGION") %>%
  mutate(
    time_slice = factor(as.character(abs(Time2)), levels = time_levels),
    covariate = factor(pretty_names[name], levels = covariate_levels),
    tree = factor(tree)
  ) %>%
  arrange(time_slice, covariate, tree)

if (nrow(coefficients) == 0) {
  stop("No coefficient rows remained after filtering the demo regression outputs.")
}

coeff_summary <- coefficients %>%
  select(
    time_slice, tree, covariate, mean, X0.025quant, X0.975quant, sig1,
    AIC, psudor2, rmse, file
  )

write.csv(
  coeff_summary,
  file = demo_here("demo_summary", "demo_coefficients_summary.csv"),
  row.names = FALSE
)

slopes <- bind_rows(lapply(slope_files, function(f) {
  d <- read_csv(f, show_col_types = FALSE)
  d$file <- basename(f)
  d
}))

slopes <- slopes %>%
  filter(col %in% names(pretty_names)) %>%
  mutate(
    time_slice = factor(as.character(yearbp), levels = time_levels),
    covariate = factor(pretty_names[col], levels = covariate_levels),
    tree = factor(tree)
  ) %>%
  arrange(time_slice, covariate, tree, x)

if (nrow(slopes) == 0) {
  stop("No slope rows remained after filtering the demo slope outputs.")
}

tree_palette <- c("1" = "#2C6E49", "2" = "#C46A2B", "3" = "#5B7FA3")
time_palette <- c("3500" = "#D1495B", "4250" = "#EDA63B", "5000" = "#2A7F9E")

format_panel_b_axis <- function(x) {
  if (length(x) > 0 && all(is.finite(x)) && max(abs(x), na.rm = TRUE) >= 1e6) {
    formatted <- paste0(
      scales::label_number(accuracy = 1, trim = TRUE)(x / 1e6),
      "M"
    )
    formatted[x == 0] <- "0"
    return(formatted)
  }

  scales::label_number(trim = TRUE)(x)
}

panel_a <- ggplot(
  coefficients,
  aes(x = mean, y = covariate, xmin = X0.025quant, xmax = X0.975quant, color = tree)
) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, color = "grey55") +
  geom_errorbar(width = 0.12, linewidth = 0.5, orientation = "y") +
  geom_point(size = 2.2) +
  facet_wrap(~ time_slice, nrow = 1) +
  scale_color_manual(values = tree_palette, name = "Tree") +
  labs(
    title = "A. Demo coefficient estimates",
    subtitle = "Three Global trees at 3500, 4250, and 5000 YBP",
    x = "Slope estimate",
    y = NULL
  ) +
  theme_cowplot(font_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold")
  )

panel_b <- ggplot(
  slopes,
  aes(x = x, y = exp(y), group = interaction(tree, time_slice), color = time_slice)
) +
  geom_line(alpha = 0.28, linewidth = 0.55) +
  facet_wrap(~ covariate, scales = "free_x", ncol = 1) +
  scale_x_continuous(labels = format_panel_b_axis) +
  scale_color_manual(values = time_palette, name = "Years BP") +
  labs(
    title = "B. Demo fitted response curves",
    subtitle = "Qualitative illustration across the 3 demo trees",
    x = "Covariate value",
    y = "Fitted clade richness"
  ) +
  theme_cowplot(font_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text = element_text(face = "bold")
  )

combined_plot <- cowplot::plot_grid(
  panel_a,
  panel_b,
  ncol = 2,
  rel_widths = c(2, 1),
  align = "h",
  axis = "tb"
)

png_path <- demo_here("demo_figures", "demo_regression_summary.png")
pdf_path <- demo_here("demo_figures", "demo_regression_summary.pdf")
session_info_path <- demo_here("demo_summary", "demo_session_info.txt")

demo_packages <- c(
  "pacman",
  "ade4", "ape", "cowplot", "data.table", "dismo", "dplyr", "ggplot2",
  "here", "INLA", "mice", "moments", "nlme", "phytools", "picante",
  "progress", "readr", "rworldmap", "scales", "sf", "sp", "spdep",
  "stringr", "tidyverse"
)

installed_demo_packages <- demo_packages[sapply(
  demo_packages,
  requireNamespace,
  quietly = TRUE
)]

package_lines <- vapply(
  sort(installed_demo_packages),
  function(pkg) {
    paste0(pkg, " ", as.character(utils::packageVersion(pkg)))
  },
  character(1)
)

ggsave(
  filename = png_path,
  plot = combined_plot,
  width = 11.5,
  height = 8,
  dpi = 300
)

ggsave(
  filename = pdf_path,
  plot = combined_plot,
  width = 11.5,
  height = 8
)

writeLines(
  c(
    "TreeSet demo environment record",
    paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "Scripts covered: code/demo/01_prepare_demo_datasets.R, code/demo/02_run_demo_regression.R, code/demo/03_plot_demo_results.R",
    "",
    "R and platform:",
    paste("R version:", R.version.string),
    paste("Platform:", R.version$platform),
    paste("Operating system:", Sys.info()[["sysname"]], Sys.info()[["release"]]),
    "",
    "Packages used across the 3 demo scripts:",
    package_lines,
    "",
    "Current plotting-session details:",
    capture.output(sessionInfo())
  ),
  con = session_info_path
)
