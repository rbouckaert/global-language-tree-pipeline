# -----------------------------------------------------------------------------#
# Supplementary MCC radial tree (gheatmap variant)
# -----------------------------------------------------------------------------#
# Purpose:
#   Build a supplementary MCC radial tree figure using ggtree::gheatmap() rings
#   so ring labels are handled as column labels (opening wedge style).
#
# This script is additive and does NOT modify the existing:
# - Fig_MCC_tree_variants_plot.R
#
# Outputs:
# - outputs/Fig_MCC_tree_variants/Fig_MCC_supp_annotations_gheatmap.(png|pdf)
# - outputs/Fig_MCC_tree_variants/Fig_MCC_supp_annotations_gheatmap_data_sources.txt
# -----------------------------------------------------------------------------#

required_plot_pkgs <- c(
  "ape", "dplyr", "ggplot2", "ggtree", "ggnewscale", "grid", "scales", "readr"
)

missing_plot_pkgs <- required_plot_pkgs[
  !vapply(required_plot_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_plot_pkgs) > 0) {
  extra_note <- if ("ggtree" %in% missing_plot_pkgs) {
    " Install ggtree via Bioconductor (e.g. BiocManager::install('ggtree'))."
  } else {
    ""
  }
  stop("Missing required plotting packages: ", paste(missing_plot_pkgs, collapse = ", "), ".", extra_note)
}
invisible(lapply(required_plot_pkgs, function(pkg) library(pkg, character.only = TRUE)))

find_data_prep_script <- function() {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- c(
    file.path(cwd, "code", "figures", "Fig_MCC_tree_variants_data_prep.R"),
    file.path(cwd, "TreeSetAnalysisScripts", "code", "figures", "Fig_MCC_tree_variants_data_prep.R"),
    file.path(cwd, "Fig_MCC_tree_variants_data_prep.R")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop(
      "Could not locate Fig_MCC_tree_variants_data_prep.R. ",
      "Run from repo root, TreeSetAnalysisScripts/, or code/figures/."
    )
  }
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

source(find_data_prep_script())

branch_colors <- c("#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#7f0000")
ed_ring_colors <- c("#deebf7", "#9ecae1", "#3182bd", "#08519c")
dr_ring_colors <- c("#e5f5e0", "#a1d99b", "#41ab5d", "#238b45", "#005a32")

base_theme <- ggplot2::theme_void() +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.key = ggplot2::element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = ggplot2::element_text(size = 7),
    legend.title = ggplot2::element_text(size = 11, face = "bold"),
    plot.caption = ggplot2::element_text(size = 8, hjust = 0),
    legend.spacing.y = grid::unit(8, "pt"),
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  )

make_branch_plot <- function(tree, branch_df, metric_col, metric_title, limits,
                             open_angle = 12, line_size = 0.22) {
  ggtree::ggtree(tree, layout = "fan", open.angle = open_angle, color = NA, size = line_size) %<+% branch_df +
    ggtree::geom_tree(ggplot2::aes(color = .data[[metric_col]]), size = line_size, lineend = "round") +
    ggplot2::scale_color_gradientn(
      colors = branch_colors,
      limits = limits,
      oob = scales::squish,
      na.value = "#BDBDBD",
      name = metric_title,
      guide = ggplot2::guide_colorbar(
        order = 1,
        title.position = "top",
        barheight = grid::unit(62, "pt"),
        barwidth = grid::unit(9, "pt")
      )
    ) +
    base_theme
}

make_ring_frame <- function(tip_df, tree_labels, value_col, ring_label) {
  if (!(value_col %in% names(tip_df))) {
    stop("Column not found in tip_df: ", value_col)
  }

  dat <- tip_df[, c("label", value_col), drop = FALSE]
  dat <- dat[match(tree_labels, dat$label), , drop = FALSE]

  if (any(is.na(dat$label))) {
    stop("Tip labels could not be aligned for ring: ", ring_label)
  }

  ring <- data.frame(value = dat[[value_col]], row.names = dat$label, check.names = FALSE)
  names(ring) <- ring_label
  ring
}

add_gheatmap_ring_cont <- function(
  p,
  ring_df,
  offset,
  width,
  legend_name,
  colors,
  limits,
  order,
  colnames_angle = 90,
  colnames_offset_y = 0.20,
  col_font_size = 2.4
) {
  x_range <- diff(range(p$data$x, na.rm = TRUE))
  if (!is.finite(x_range) || x_range <= 0) {
    stop("Cannot determine plot x-range for gheatmap width scaling.")
  }
  width_scaled <- width / x_range

  p <- ggtree::gheatmap(
    p,
    ring_df,
    offset = offset,
    width = width_scaled,
    color = NA,
    colnames = TRUE,
    colnames_position = "top",
    colnames_angle = colnames_angle,
    colnames_offset_y = colnames_offset_y,
    font.size = col_font_size
  )

  p + ggplot2::scale_fill_gradientn(
    colors = colors,
    limits = limits,
    oob = scales::squish,
    na.value = "#F8F8F8",
    name = legend_name,
    guide = ggplot2::guide_colorbar(
      order = order,
      title.position = "top",
      barheight = grid::unit(36, "pt"),
      barwidth = grid::unit(9, "pt")
    )
  )
}

add_gheatmap_ring_disc <- function(
  p,
  ring_df,
  offset,
  width,
  legend_name,
  values,
  order,
  colnames_angle = 90,
  colnames_offset_y = 0.20,
  col_font_size = 2.4,
  ncol = 1,
  drop = TRUE
) {
  x_range <- diff(range(p$data$x, na.rm = TRUE))
  if (!is.finite(x_range) || x_range <= 0) {
    stop("Cannot determine plot x-range for gheatmap width scaling.")
  }
  width_scaled <- width / x_range

  p <- ggtree::gheatmap(
    p,
    ring_df,
    offset = offset,
    width = width_scaled,
    color = NA,
    colnames = TRUE,
    colnames_position = "top",
    colnames_angle = colnames_angle,
    colnames_offset_y = colnames_offset_y,
    font.size = col_font_size
  )

  p + ggplot2::scale_fill_manual(
    values = values,
    na.value = "#F8F8F8",
    drop = drop,
    name = legend_name,
    guide = ggplot2::guide_legend(
      order = order,
      ncol = ncol,
      override.aes = list(color = NA)
    )
  )
}

save_plot_pair <- function(plot_obj, out_dir, stem, width = 14.0, height = 11.8, scale = 1) {
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"

  ggplot2::ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    scale = scale
  )
  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white",
    scale = scale
  )
}

build_supp_gheatmap_plot <- function(prepped) {
  tree <- prepped$tree
  tip_df <- prepped$tip_df

  p <- make_branch_plot(
    tree = tree,
    branch_df = prepped$branch_ed,
    metric_col = "branch_ED_clip",
    metric_title = "Branch ED",
    limits = prepped$ed_branch_clip$limits,
    open_angle = 12,
    line_size = 0.22
  )

  # Ring geometry in tree-x units; gheatmap uses absolute x coordinates.
  n_tip <- ape::Ntip(tree)
  tip_depth <- max(ape::node.depth.edgelength(tree)[seq_len(n_tip)], na.rm = TRUE)
  ring_offset <- tip_depth * 0.002
  ring_width <- tip_depth * 0.040
  ring_gap <- tip_depth * 0.000
  col_offset_y <- ring_width * 0.35

  ring_order <- 3
  selected_cont_ordered <- c("popd_log", "dist_city_log", "friction_log", "area_log")
  selected_cont_ordered <- selected_cont_ordered[selected_cont_ordered %in% prepped$selected_cont]

  cont_specs <- list(
    popd_log = list(
      ring_label = "PopD",
      legend_name = "Population density (log)",
      colors = c("#f2f0f7", "#9e9ac8", "#54278f")
    ),
    dist_city_log = list(
      ring_label = "City",
      legend_name = "Distance to city (log)",
      colors = c("#f7fbff", "#6baed6", "#08306b")
    ),
    friction_log = list(
      ring_label = "Fric",
      legend_name = "Friction (log)",
      colors = c("#fff7bc", "#d9a441", "#6b4226")
    ),
    area_log = list(
      ring_label = "Area",
      legend_name = "Area (log)",
      colors = c("#f7fcf0", "#74c476", "#00441b")
    )
  )

  # 1) Tip ED ring
  ring_df <- make_ring_frame(tip_df, tree$tip.label, "ED_tip_clip", "ED")
  p <- add_gheatmap_ring_cont(
    p = p,
    ring_df = ring_df,
    offset = ring_offset,
    width = ring_width,
    legend_name = "Tip ED",
    colors = ed_ring_colors,
    limits = prepped$ed_tip_clip$limits,
    order = ring_order,
    colnames_offset_y = col_offset_y
  )
  ring_order <- ring_order + 1
  ring_offset <- ring_offset + ring_width + ring_gap

  # 2) Optional DR ring
  if (isTRUE(prepped$show_dr_ring)) {
    p <- p + ggnewscale::new_scale_fill()
    ring_df <- make_ring_frame(tip_df, tree$tip.label, "DR_tip_clip", "DR")
    p <- add_gheatmap_ring_cont(
      p = p,
      ring_df = ring_df,
      offset = ring_offset,
      width = ring_width,
      legend_name = "Tip DR (1/ED)",
      colors = dr_ring_colors,
      limits = prepped$dr_tip_clip$limits,
      order = ring_order,
      colnames_offset_y = col_offset_y
    )
    ring_order <- ring_order + 1
    ring_offset <- ring_offset + ring_width + ring_gap
  }

  # 3) Continuous predictor rings (fixed narrative order)
  for (pred_name in selected_cont_ordered) {
    pred_clip <- clip_quantiles(tip_df[[pred_name]], probs = c(0.01, 0.99))
    clip_col <- paste0(pred_name, "_clip")
    tip_df[[clip_col]] <- pred_clip$values

    p <- p + ggnewscale::new_scale_fill()
    ring_df <- make_ring_frame(tip_df, tree$tip.label, clip_col, cont_specs[[pred_name]]$ring_label)
    p <- add_gheatmap_ring_cont(
      p = p,
      ring_df = ring_df,
      offset = ring_offset,
      width = ring_width,
      legend_name = cont_specs[[pred_name]]$legend_name,
      colors = cont_specs[[pred_name]]$colors,
      limits = pred_clip$limits,
      order = ring_order,
      colnames_offset_y = col_offset_y
    )

    ring_order <- ring_order + 1
    ring_offset <- ring_offset + ring_width + ring_gap
  }

  # 4) Optional island ring
  if (isTRUE(prepped$show_island)) {
    p <- p + ggnewscale::new_scale_fill()
    ring_df <- make_ring_frame(tip_df, tree$tip.label, "island_bin", "Island")
    p <- add_gheatmap_ring_disc(
      p = p,
      ring_df = ring_df,
      offset = ring_offset,
      width = ring_width,
      legend_name = "Island",
      values = c(Mainland = "#d9d9d9", Island = "#1f78b4"),
      order = ring_order,
      colnames_offset_y = col_offset_y,
      ncol = 1,
      drop = FALSE
    )

    ring_order <- ring_order + 1
    ring_offset <- ring_offset + ring_width + ring_gap
  }

  # 5) Family ring (outermost)
  p <- p + ggnewscale::new_scale_fill()
  ring_df <- make_ring_frame(tip_df, tree$tip.label, "top_family", "Family")
  p <- add_gheatmap_ring_disc(
    p = p,
    ring_df = ring_df,
    offset = ring_offset,
    width = ring_width,
    legend_name = "Top 27 families",
    values = prepped$family_palette,
    order = ring_order,
    colnames_offset_y = col_offset_y,
    ncol = 3,
    drop = FALSE
  )

  supp_tracks <- c("tip ED")
  if (isTRUE(prepped$show_dr_ring)) supp_tracks <- c(supp_tracks, "tip DR")
  if (length(selected_cont_ordered) > 0) {
    supp_tracks <- c(supp_tracks, gsub("_log", "", selected_cont_ordered, fixed = TRUE))
  }
  if (isTRUE(prepped$show_island)) supp_tracks <- c(supp_tracks, "island")
  supp_tracks <- c(supp_tracks, "family")

  p <- p + ggplot2::labs(
    caption = paste0("Tracks shown: ", paste(supp_tracks, collapse = ", "))
  )

  list(plot = p, tip_df = tip_df, supp_tracks = supp_tracks)
}

write_supp_gheatmap_outputs <- function(prepped, built) {
  out_dir <- prepped$paths$out_dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  save_plot_pair(
    plot_obj = built$plot,
    out_dir = out_dir,
    stem = "Fig_MCC_supp_annotations_gheatmap",
    width = 14.0,
    height = 11.8
  )

  source_lines <- c(
    paste0("MCC tree: ", prepped$paths$tree_file),
    paste0("Top27 families: ", prepped$paths$top27_file),
    paste0("ED source used: ", prepped$metric_info$ed_source),
    paste0("DR source used: ", prepped$metric_info$dr_source),
    paste0("Figure file stem: Fig_MCC_supp_annotations_gheatmap"),
    paste0(
      "Predictor source(s): ",
      if (length(prepped$pred_info$sources) > 0) paste(prepped$pred_info$sources, collapse = "; ") else "none detected"
    ),
    paste0(
      "Predictor Time3 filter: Time3 == ", prepped$target_time3,
      " (applied to sources containing Time3)"
    ),
    paste0(
      "Predictor inclusion threshold: ",
      round(100 * prepped$predictor_coverage_threshold), "% non-missing tips"
    ),
    paste0(
      "Selected continuous predictors: ",
      if (length(prepped$selected_cont) > 0) paste(prepped$selected_cont, collapse = ", ") else "none"
    ),
    paste0("Island strip included: ", prepped$show_island),
    paste0("DR strip included (only when no continuous predictors selected): ", prepped$show_dr_ring),
    paste0("Tracks shown in plot: ", paste(built$supp_tracks, collapse = ", ")),
    "Ring engine: ggtree::gheatmap with ggnewscale::new_scale_fill"
  )

  writeLines(
    source_lines,
    con = file.path(out_dir, "Fig_MCC_supp_annotations_gheatmap_data_sources.txt")
  )

  message("Saved gheatmap supplementary figure outputs in: ", out_dir)
}

has_valid_prepped <- function(x) {
  if (!is.list(x)) return(FALSE)

  required_fields <- c(
    "paths", "tree", "tip_df", "branch_ed", "ed_branch_clip", "ed_tip_clip",
    "dr_tip_clip", "family_palette", "family_labels", "pred_info", "selected_cont",
    "show_island", "show_dr_ring", "metric_info", "target_time3",
    "predictor_coverage_threshold"
  )
  all(required_fields %in% names(x))
}

reuse_prepped <- exists("prepped", envir = .GlobalEnv, inherits = FALSE) &&
  has_valid_prepped(get("prepped", envir = .GlobalEnv, inherits = FALSE))

if (reuse_prepped) {
  message("Using existing `prepped` from .GlobalEnv; skipping data loading.")
  prepped <- get("prepped", envir = .GlobalEnv, inherits = FALSE)
} else {
  prepped <- prepare_mcc_tree_variant_data(
    base_dir = resolve_base_dir(),
    target_time3 = TARGET_TIME3,
    predictor_coverage_threshold = 0.80
  )
  assign("prepped", prepped, envir = .GlobalEnv)
  message("Created and cached `prepped` in .GlobalEnv.")
}

built <- build_supp_gheatmap_plot(prepped)
write_supp_gheatmap_outputs(prepped, built)
