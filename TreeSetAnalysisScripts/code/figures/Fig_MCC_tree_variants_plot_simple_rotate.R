# -----------------------------------------------------------------------------#
# Figure script: MCC radial tree variants (main + supplementary)
# -----------------------------------------------------------------------------#
# Purpose:
#   Build and save the MCC tree figures using preprocessed data from
#   Fig_MCC_tree_variants_data_prep.R.
#
# Run from:
# - TreeSetAnalysisScripts/ (recommended)
# - repository root
# -----------------------------------------------------------------------------#

required_plot_pkgs <- c(
  "ggplot2", "ggtree", "ggtreeExtra", "ggnewscale", "grid", "scales", "readr", "dplyr"
)

missing_plot_pkgs <- required_plot_pkgs[
  !vapply(required_plot_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_plot_pkgs) > 0) {
  extra_note <- if ("ggtree" %in% missing_plot_pkgs || "ggtreeExtra" %in% missing_plot_pkgs) {
    " Install ggtree/ggtreeExtra via Bioconductor (e.g. BiocManager::install(c('ggtree','ggtreeExtra')))."
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
    plot.title = ggplot2::element_text(size = 13, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 10),
    plot.caption = ggplot2::element_text(size = 8, hjust = 0),
    legend.spacing.y = grid::unit(8, "pt"),
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  )

make_branch_plot <- function(tree, branch_df, metric_col, metric_title, limits,
                             open_angle = 0, line_size = 0.18) {
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

add_family_ring <- function(p, tip_df, family_palette, family_labels,
                            offset = 0.015, pwidth = 0.045, width = 1,
                            add_new_scale = FALSE) {
  if (isTRUE(add_new_scale)) {
    p <- p + ggnewscale::new_scale_fill()
  }

  p +
    ggtreeExtra::geom_fruit(
      data = tip_df,
      geom = geom_tile,
      mapping = ggplot2::aes(y = label, x = 1, fill = top_family),
      offset = offset,
      pwidth = pwidth,
      width = width,
      color = NA
    ) +
    ggplot2::scale_fill_manual(
      values = family_palette,
      labels = family_labels,
      drop = FALSE,
      name = "Family",
      guide = ggplot2::guide_legend(
        order = 2,
        ncol = 7,
        byrow = TRUE,
        position = "bottom",
        keywidth = grid::unit(10, "pt"),
        keyheight = grid::unit(10, "pt"),
        theme = ggplot2::theme(legend.spacing.y = grid::unit(2, "pt")),
        override.aes = list(color = NA)
      )
    )
}

add_major_family_labels <- function(p, family_nodes, offset = 0.022, offset_text = 0.012, fontsize = 2.8) {
  if (nrow(family_nodes) == 0) return(p)

  p + ggtree::geom_cladelab(
    data = family_nodes,
    mapping = ggplot2::aes(node = node, label = label),
    align = TRUE,
    offset = offset,
    offset.text = offset_text,
    barsize = 0.2,
    fontsize = fontsize,
    angle = "auto",
    horizontal = FALSE,
    color = "#2F2F2F"
  )
}

add_major_family_arcs <- function(
  p,
  family_nodes,
  offset = 40,
  offset_text = 10,
  barsize = 0.2,
  color = "#2F2F2F"
) {
  if (nrow(family_nodes) == 0) return(p)
  arc_df <- family_nodes
  arc_df$label <- ""

  p + ggtree::geom_cladelab(
    data = arc_df,
    mapping = ggplot2::aes(node = node, label = label),
    align = TRUE,
    offset = offset,
    offset.text = offset_text,
    barsize = barsize,
    fontsize = 0,
    angle = "auto",
    horizontal = FALSE,
    color = color
  )
}

add_family_arc_labels <- function(p, label_df, offset = 0.07, pwidth = 0.08, size = 2.8) {
  if (nrow(label_df) == 0) return(p)

  p + ggtreeExtra::geom_fruit(
    data = label_df,
    geom = geom_text,
    mapping = ggplot2::aes(
      y = label,
      x = 1,
      label = ring_label,
      angle = text_angle,
      hjust = text_hjust
    ),
    offset = offset,
    pwidth = pwidth,
    color = "#2F2F2F",
    size = size,
    vjust = 0.5,
    inherit.aes = FALSE
  )
}

compute_tangent_family_label_angles <- function(p, label_df) {
  if (is.null(label_df) || nrow(label_df) == 0) return(label_df)
  if (!all(c("label", "ring_label") %in% names(label_df))) return(label_df)

  tip_pos <- p$data[p$data$isTip %in% TRUE, c("label", "y"), drop = FALSE]
  names(tip_pos)[names(tip_pos) == "y"] <- "tip_y"
  if (nrow(tip_pos) == 0) return(label_df[0, , drop = FALSE])

  out <- dplyr::left_join(label_df, tip_pos, by = "label")
  out <- out[is.finite(out$tip_y), , drop = FALSE]
  if (nrow(out) == 0) return(out)

  pb <- ggplot2::ggplot_build(p)
  theta_range <- pb$layout$panel_params[[1]]$theta.range
  if (length(theta_range) != 2 || !all(is.finite(theta_range))) return(out)

  span <- diff(theta_range)
  if (!is.finite(span) || span <= 0) return(out)

  start <- pb$plot$coordinates$start
  direction <- pb$plot$coordinates$direction

  theta <- (start + direction * 2 * pi * ((out$tip_y - theta_range[1]) / span)) %% (2 * pi)
  text_angle <- -theta * 180 / pi

  text_hjust <- ifelse(text_angle < -90 | text_angle > 90, 1, 0)
  text_angle <- ifelse(
    text_angle < -90,
    text_angle + 180,
    ifelse(text_angle > 90, text_angle - 180, text_angle)
  )

  out$text_angle <- text_angle
  out$text_hjust <- text_hjust
  out$tip_y <- NULL
  out
}

add_opening_ring_labels <- function(
  p,
  tree,
  ring_label_df,
  open_angle = 12,
  gap_position = c("right", "top", "left", "bottom"),
  size = 1.8,
  color = "#4A4A4A"
) {
  if (nrow(ring_label_df) == 0) return(p)
  gap_position <- match.arg(gap_position)

  pb <- ggplot2::ggplot_build(p)

  is_tile <- vapply(pb$data, function(d) {
    all(c("xmin", "xmax", "ymin", "ymax") %in% names(d)) && nrow(d) > 500
  }, logical(1))
  tile_idx <- which(is_tile)

  n <- min(nrow(ring_label_df), length(tile_idx))
  if (n == 0) return(p)
  ring_label_df <- ring_label_df[seq_len(n), , drop = FALSE]

  ring_x <- vapply(tile_idx[seq_len(n)], function(i) {
    median(pb$data[[i]]$x, na.rm = TRUE)
  }, numeric(1))

  # Calculate gap center for y
  n_tip <- ape::Ntip(tree)
  gap_size_y <- 1 + n_tip * open_angle / (360 - open_angle)
  y_center <- n_tip + gap_size_y / 2

  text_angle <- switch(
    gap_position,
    right = 270,  # Tangent orientation at 3 o'clock
    top = 0,      # Tangent orientation at 12 o'clock
    left = 90,    # Tangent orientation at 9 o'clock
    bottom = 180  # Tangent orientation at 6 o'clock
  )

  label_df <- data.frame(
    x_label = ring_x,
    y_label = y_center,
    ring_label = ring_label_df$ring_label,
    stringsAsFactors = FALSE
  )

  p + ggplot2::geom_text(
    data = label_df,
    mapping = ggplot2::aes(x = x_label, y = y_label, label = ring_label),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    angle = text_angle,
    size = size,
    color = color
  )
}

add_continuous_ring <- function(p, tip_df, value_col, legend_name, colors, limits,
                                offset, pwidth, order, width = 1) {
  p +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data = tip_df,
      geom = geom_tile,
      mapping = ggplot2::aes(y = label, x = 1, fill = .data[[value_col]]),
      offset = offset,
      pwidth = pwidth,
      width = width,
      color = NA
    ) +
    ggplot2::scale_fill_gradientn(
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

add_binary_ring <- function(p, tip_df, value_col, legend_name, values,
                            offset, pwidth, order, width = 1) {
  p +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data = tip_df,
      geom = geom_tile,
      mapping = ggplot2::aes(y = label, x = 1, fill = .data[[value_col]]),
      offset = offset,
      pwidth = pwidth,
      width = width,
      color = NA
    ) +
    ggplot2::scale_fill_manual(
      values = values,
      na.value = "#EFEFEF",
      name = legend_name,
      guide = ggplot2::guide_legend(
        order = order,
        override.aes = list(color = NA)
      )
    )
}

save_plot_pair <- function(plot_obj, out_dir, stem, width = 12, height = 11, scale = 1, filename_suffix = "") {
  png_file <- file.path(out_dir, paste0(stem, filename_suffix, ".png"))
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
  # ggplot2::ggsave(
  #   filename = pdf_file,
  #   plot = plot_obj,
  #   width = width,
  #   height = height,
  #   units = "in",
  #   device = pdf_device,
  #   bg = "white",
  #   scale = scale
  # )
}

default_plot_tuning <- function() {
  list(
    main = list(
      open_angle = 0,
      line_size = 0.19,
      family_ring_offset = 0.014,
      family_ring_pwidth = 0.090,
      family_ring_width = 3,
      family_label_offset = 8.0,
      family_label_offset_text = 7.0,
      family_label_fontsize = 3.2
    ),
    supp = list(
      gap_position = "top",
      open_angle = 30,
      line_size = 0.22,
      ring_order_start = 3,
      ring_offset_start = 0,
      ring_pwidth = 0.040,
      ring_tile_width = 7,
      family_ring_offset = NULL,
      family_ring_pwidth = NULL,
      family_ring_width = 7,
      rotation_angle = 102,
      family_label_offset = 40,
      family_label_offset_text = 10,
      family_label_fontsize = 3.4,
      austronesian_angle = -72,
      opening_ring_label_size = 3
    )
  )
}

build_mcc_tree_variant_plots <- function(prepped, plot_tuning = list()) {
  tuning <- utils::modifyList(default_plot_tuning(), plot_tuning)
  main_tuning <- tuning$main
  supp_tuning <- tuning$supp

  tree <- prepped$tree
  tip_df <- prepped$tip_df

  p_main <- make_branch_plot(
    tree = tree,
    branch_df = prepped$branch_ed,
    metric_col = "branch_ED_clip",
    metric_title = "Branch ED",
    limits = prepped$ed_branch_clip$limits,
    open_angle = main_tuning$open_angle,
    line_size = main_tuning$line_size
  )

  p_main <- add_family_ring(
    p = p_main,
    tip_df = tip_df,
    family_palette = prepped$family_palette,
    family_labels = prepped$family_labels,
    offset = main_tuning$family_ring_offset,
    pwidth = main_tuning$family_ring_pwidth,
    width = main_tuning$family_ring_width,
    add_new_scale = FALSE
  )

  p_main <- add_major_family_labels(
    p = p_main,
    family_nodes = prepped$family_nodes,
    offset = main_tuning$family_label_offset,
    offset_text = main_tuning$family_label_offset_text,
    fontsize = main_tuning$family_label_fontsize
  )

  # p_main <- p_main +
  #   ggplot2::ggtitle("MCC radial tree: branch ED and top 27 families") +
  #   ggplot2::labs(
  #     caption = paste0(
  #       "ED source: ", prepped$metric_info$ed_source,
  #       "; branch colour clipped to 1st-99th percentiles."
  #     )
  #   )

  # Supplementry ------------------
  supp_gap_position <- supp_tuning$gap_position
  supp_open_angle <- supp_tuning$open_angle

  p_supp <- make_branch_plot(
    tree = tree,
    branch_df = prepped$branch_ed,
    metric_col = "branch_ED_clip",
    metric_title = "Branch ED",
    limits = prepped$ed_branch_clip$limits,
    open_angle = supp_open_angle,
    line_size = supp_tuning$line_size
  )

  ring_order <- supp_tuning$ring_order_start
  offset <- supp_tuning$ring_offset_start
  ring_pwidth_cont <- supp_tuning$ring_pwidth
  ring_tile_width <- supp_tuning$ring_tile_width
  ring_label_df <- data.frame(
    ring_label = character(0),
    ring_offset = numeric(0),
    ring_pwidth = numeric(0),
    stringsAsFactors = FALSE
  )

  append_ring_label <- function(df, ring_label, ring_offset, ring_pwidth) {
    rbind(
      df,
      data.frame(
        ring_label = ring_label,
        ring_offset = ring_offset,
        ring_pwidth = ring_pwidth,
        stringsAsFactors = FALSE
      )
    )
  }

  p_supp <- add_continuous_ring(
    p = p_supp,
    tip_df = tip_df,
    value_col = "ED_tip_clip",
    legend_name = "Tip ED",
      colors = ed_ring_colors,
      limits = prepped$ed_tip_clip$limits,
      offset = 0,
      pwidth = ring_pwidth_cont,
      order = ring_order,
      width = ring_tile_width
    )
  ring_label_df <- append_ring_label(ring_label_df, "ED", offset, ring_pwidth_cont)
  pwidth_curr <- ring_pwidth_cont
  ring_order <- ring_order + 1
  offset <- offset

  if (isTRUE(prepped$show_dr_ring)) {
    p_supp <- add_continuous_ring(
      p = p_supp,
      tip_df = tip_df,
      value_col = "DR_tip_clip",
      legend_name = "Tip DR (1/ED)",
      colors = dr_ring_colors,
      limits = prepped$dr_tip_clip$limits,
      offset = offset,
      pwidth = ring_pwidth_cont,
      order = ring_order,
      width = ring_tile_width
    )
    ring_label_df <- append_ring_label(ring_label_df, "DR", offset, ring_pwidth_cont)
    pwidth_curr <- ring_pwidth_cont
    ring_order <- ring_order + 1
    offset <- offset
  }

  cont_specs <- list(
    area_log = list(name = "Area (log)", colors = c("#f7fcf0", "#74c476", "#00441b")),
    popd_log = list(name = "Population density (log)", colors = c("#f2f0f7", "#9e9ac8", "#54278f")),
    dist_city_log = list(name = "Distance to city (log)", colors = c("#f7fbff", "#6baed6", "#08306b")),
    friction_log = list(name = "Friction (log)", colors = c("#fff7bc", "#d9a441", "#6b4226"))
  )
  cont_short_labels <- c(
    popd_log = "PopD",
    dist_city_log = "City",
    friction_log = "Fric",
    area_log = "Area"
  )
  selected_cont_ordered <- c("popd_log", "dist_city_log", "friction_log", "area_log")
  selected_cont_ordered <- selected_cont_ordered[selected_cont_ordered %in% prepped$selected_cont]

  for (pred_name in selected_cont_ordered) {
    pred_clip <- clip_quantiles(tip_df[[pred_name]], probs = c(0.01, 0.99))
    tip_df[[paste0(pred_name, "_clip")]] <- pred_clip$values

    p_supp <- add_continuous_ring(
      p = p_supp,
      tip_df = tip_df,
      value_col = paste0(pred_name, "_clip"),
      legend_name = cont_specs[[pred_name]]$name,
      colors = cont_specs[[pred_name]]$colors,
      limits = pred_clip$limits,
      offset = offset,
      pwidth = ring_pwidth_cont,
      order = ring_order,
      width = ring_tile_width
    )
    ring_label_df <- append_ring_label(ring_label_df, cont_short_labels[[pred_name]], offset, ring_pwidth_cont)
    pwidth_curr <- ring_pwidth_cont
    ring_order <- ring_order + 1
    offset <- offset
  }

  if (isTRUE(prepped$show_island)) {
    p_supp <- add_binary_ring(
      p = p_supp,
      tip_df = tip_df,
      value_col = "island_bin",
      legend_name = "Island",
      values = c(Mainland = "#d9d9d9", Island = "#1f78b4"),
      offset = offset,
      pwidth = ring_pwidth_cont,
      order = ring_order,
      width = ring_tile_width
    )
    ring_label_df <- append_ring_label(ring_label_df, "Island", offset, ring_pwidth_cont)
    pwidth_curr <- ring_pwidth_cont
    ring_order <- ring_order + 1
    offset <- offset
  }

  # Keep family as the outermost annotation ring in the supplementary panel.
  family_pwidth <- if (is.null(supp_tuning$family_ring_pwidth)) ring_pwidth_cont else supp_tuning$family_ring_pwidth
  family_offset <- if (is.null(supp_tuning$family_ring_offset)) offset else supp_tuning$family_ring_offset
  p_supp <- add_family_ring(
    p = p_supp,
    tip_df = tip_df,
    family_palette = prepped$family_palette,
    family_labels = prepped$family_labels,
    offset = family_offset,
    pwidth = family_pwidth,
    width = supp_tuning$family_ring_width,
    add_new_scale = TRUE
  )
  # Simple explicit rotation knob for the supplementary fan plot.
  # Positive values rotate clockwise; negative values rotate counterclockwise.
  supp_rotation_angle <- supp_tuning$rotation_angle
  if (abs(supp_rotation_angle) > 1e-8) {
    p_supp <- ggtree::rotate_tree(p_supp, angle = supp_rotation_angle)
  }

  # Keep auto cladelabel placement (works for nearly all families), and patch only
  # the one label that can flip upside-down after rotation.
  supp_family_nodes <- prepped$family_nodes
  supp_family_nodes$label[supp_family_nodes$label == "Austronesian"] <- ""

  p_supp <- add_major_family_labels(
    p = p_supp,
    family_nodes = supp_family_nodes,
    offset = supp_tuning$family_label_offset,
    offset_text = supp_tuning$family_label_offset_text,
    fontsize = supp_tuning$family_label_fontsize
  )

  austro_row <- prepped$family_nodes[prepped$family_nodes$label == "Austronesian", , drop = FALSE]
  if (nrow(austro_row) > 0) {
    p_supp <- p_supp + ggtree::geom_cladelab(
      data = austro_row,
      mapping = ggplot2::aes(node = node, label = label),
      align = TRUE,
      offset = supp_tuning$family_label_offset,
      offset.text = supp_tuning$family_label_offset_text,
      barsize = 0.2,
      fontsize = supp_tuning$family_label_fontsize,
      angle = supp_tuning$austronesian_angle,
      horizontal = FALSE,
      color = "#2F2F2F"
    )
  }

  p_supp <- add_opening_ring_labels(
    p = p_supp,
    tree = tree,
    ring_label_df = ring_label_df,
    open_angle = supp_open_angle,
    gap_position = supp_gap_position,
    size = supp_tuning$opening_ring_label_size
  )

  supp_tracks <- c("tip ED")
  if (isTRUE(prepped$show_dr_ring)) supp_tracks <- c(supp_tracks, "tip DR")
  if (length(selected_cont_ordered) > 0) {
    supp_tracks <- c(supp_tracks, gsub("_log", "", selected_cont_ordered, fixed = TRUE))
  }
  if (isTRUE(prepped$show_island)) supp_tracks <- c(supp_tracks, "island")
  supp_tracks <- c(supp_tracks, "family")

  # p_supp <- p_supp +
  #   ggplot2::labs(
  #     caption = paste0(
  #       "Tracks shown: ", paste(supp_tracks, collapse = ", ")
  #     )
  #   )

  list(
    p_main = p_main,
    p_supp = p_supp,
    supp_tracks = supp_tracks,
    tip_df = tip_df,
    supp_gap_position = supp_gap_position,
    supp_rotation_angle = supp_rotation_angle,
    supp_open_angle = supp_open_angle,
    plot_tuning = tuning
  )
}

write_mcc_tree_variant_outputs <- function(prepped, built_plots) {
  out_dir <- prepped$paths$out_dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  angle_suffix <- paste0("_", as.integer(round(built_plots$supp_rotation_angle)), "degree")

  save_plot_pair(
    built_plots$p_main,
    out_dir,
    "Fig_MCC_main_topology_ED",
    width = 11.6,
    height = 10.6,
    filename_suffix = angle_suffix
  )
  save_plot_pair(
    built_plots$p_supp,
    out_dir,
    "Fig_MCC_supp_annotations_simple_rotate",
    width = 14.0,
    height = 11.8,
    filename_suffix = angle_suffix
  )

  tip_df <- built_plots$tip_df
  tip_export <- tip_df %>%
    dplyr::mutate(
      family_label = pretty_family_name(as.character(top_family)),
      family_is_top27 = top_family != "Other"
    ) %>%
    dplyr::select(
      label, top_family, family_label, family_is_top27,
      ED_tip, DR_tip, ED_tip_clip, DR_tip_clip,
      area_raw, popd_raw, dist_city_raw, friction_raw, island_raw,
      area_log, popd_log, dist_city_log, friction_log, island_bin
    )

  #readr::write_csv(tip_export, file.path(out_dir, "Fig_MCC_tree_variants_tip_metadata.csv"))

  if (nrow(prepped$family_nodes) > 0) {
    #readr::write_csv(prepped$family_nodes, file.path(out_dir, "Fig_MCC_tree_variants_family_labels.csv"))
  }

  coverage_lines <- c(
    sprintf("area_log coverage: %.3f", coverage_prop(tip_df$area_log)),
    sprintf("popd_log coverage: %.3f", coverage_prop(tip_df$popd_log)),
    sprintf("dist_city_log coverage: %.3f", coverage_prop(tip_df$dist_city_log)),
    sprintf("friction_log coverage: %.3f", coverage_prop(tip_df$friction_log)),
    sprintf("island_bin coverage: %.3f", coverage_prop(tip_df$island_bin))
  )

  source_lines <- c(
    paste0("MCC tree: ", prepped$paths$tree_file),
    paste0("Top27 families: ", prepped$paths$top27_file),
    paste0("ED source used: ", prepped$metric_info$ed_source),
    paste0("DR source used: ", prepped$metric_info$dr_source),
    paste0("Main figure file stem: Fig_MCC_main_topology_ED", angle_suffix),
    paste0("Supplementary figure file stem: Fig_MCC_supp_annotations_simple_rotate", angle_suffix),
    paste0("ED branch clipping limits (1%-99%): ", paste(round(prepped$ed_branch_clip$limits, 4), collapse = " to ")),
    paste0("ED tip clipping limits (1%-99%): ", paste(round(prepped$ed_tip_clip$limits, 4), collapse = " to ")),
    paste0("DR tip clipping limits (1%-99%): ", paste(round(prepped$dr_tip_clip$limits, 6), collapse = " to ")),
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
    paste0("Family ring set as outermost in supplementary figure: TRUE"),
    paste0("Supplementary opening angle: ", built_plots$supp_open_angle, " degrees"),
    paste0(
      "Supplementary opening position: ",
      built_plots$supp_gap_position,
      " (rotation angle ",
      built_plots$supp_rotation_angle,
      " degrees)"
    ),
    coverage_lines
  )

  writeLines(source_lines, con = file.path(out_dir, "Fig_MCC_tree_variants_data_sources.txt"))
  message("Saved figure outputs and metadata in: ", out_dir)
}

has_valid_prepped <- function(x) {
  if (!is.list(x)) return(FALSE)

  required_fields <- c(
    "base_dir", "target_time3", "predictor_coverage_threshold",
    "paths", "tree", "tip_df", "branch_ed", "branch_dr",
    "ed_branch_clip", "dr_branch_clip", "ed_tip_clip", "dr_tip_clip",
    "family_nodes", "family_ring_labels", "family_palette", "family_labels",
    "pred_info", "pred_coverage", "selected_cont", "show_island",
    "show_dr_ring", "metric_info"
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

# All major plot tuning knobs live here so you can adjust in one place.
plot_tuning <- list(
  main = list(
    open_angle = 0,
    line_size = 0.19,
    family_ring_offset = 0.014,
    family_ring_pwidth = 0.090,
    family_ring_width = 3,
    family_label_offset = 8.0,
    family_label_offset_text = 7.0,
    family_label_fontsize = 3.2
  ),
  supp = list(
    gap_position = "top",
    open_angle = 30,
    line_size = 0.22,
    ring_order_start = 3,
    ring_offset_start = 0.03,
    ring_pwidth = 0.040,
    ring_tile_width = 7,
    family_ring_offset = 0.1,
    family_ring_pwidth = NULL,
    family_ring_width = 7,
    rotation_angle = 102,
    family_label_offset = 65,
    family_label_offset_text = 10,
    family_label_fontsize = 3.4,
    austronesian_angle = -50,
    opening_ring_label_size = 3.5
  )
)

built_plots <- build_mcc_tree_variant_plots(prepped, plot_tuning = plot_tuning)
write_mcc_tree_variant_outputs(prepped, built_plots)
