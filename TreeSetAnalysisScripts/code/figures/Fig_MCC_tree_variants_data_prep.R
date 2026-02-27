# -----------------------------------------------------------------------------#
# Data prep for MCC radial tree variants
# -----------------------------------------------------------------------------#
# Purpose:
#   Centralize all data loading, tree cleaning, metric computation, and
#   predictor merging used by MCC tree variant figures.
#
# Exposes:
# - ensure_mcc_data_packages()
# - resolve_base_dir()
# - prepare_mcc_tree_variant_data()
# -----------------------------------------------------------------------------#

TARGET_TIME3 <- -4250

ensure_mcc_data_packages <- function() {
  required_pkgs <- c("ape", "dplyr", "picante", "phytools", "readr", "tidyr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop("Missing required packages for data prep: ", paste(missing_pkgs, collapse = ", "), ".")
  }

  invisible(lapply(required_pkgs, function(pkg) library(pkg, character.only = TRUE)))
}

resolve_base_dir <- function() {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

  if (basename(cwd) == "TreeSetAnalysisScripts") {
    return(cwd)
  }

  candidate <- file.path(cwd, "TreeSetAnalysisScripts")
  if (dir.exists(candidate)) {
    return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
  }

  stop("Could not locate TreeSetAnalysisScripts. Run from repo root or TreeSetAnalysisScripts/.")
}

normalize_glottocode <- function(x) {
  x <- trimws(as.character(x))
  x[x == "osse1243"] <- "iron1242"
  x
}

pretty_family_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("Bantu (Atlantic-Congo)", "Bantu-Atlantic-Congo", x, fixed = TRUE)
  x <- gsub("Central-Sudanic", "Central Sudanic", x, fixed = TRUE)
  x <- gsub("North-America", "North America", x, fixed = TRUE)
  x <- gsub("South-America", "South America", x, fixed = TRUE)
  x <- gsub("Nuclear-Trans-New-Guinea", "Nuclear Trans New Guinea", x, fixed = TRUE)
  x <- gsub("Nuclear-Torricelli", "Nuclear Torricelli", x, fixed = TRUE)
  x
}

first_present_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

first_non_na <- function(x) {
  idx <- which(!is.na(x))
  if (length(idx) == 0) NA_real_ else x[idx[1]]
}

pull_numeric_col <- function(df, candidates) {
  col_nm <- first_present_col(df, candidates)
  if (is.na(col_nm)) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[col_nm]]))
}

coverage_prop <- function(x) {
  if (length(x) == 0) return(0)
  mean(!is.na(x))
}

clip_quantiles <- function(x, probs = c(0.01, 0.99)) {
  x <- as.numeric(x)
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0) {
    return(list(values = x, limits = c(NA_real_, NA_real_)))
  }

  q <- stats::quantile(finite_x, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
  if (!all(is.finite(q)) || q[1] >= q[2]) {
    q <- range(finite_x, na.rm = TRUE)
  }

  clipped <- pmax(pmin(x, q[2]), q[1])
  list(values = clipped, limits = q)
}

read_mcc_tree_safe <- function(tree_file) {
  tree_try <- tryCatch(ape::read.nexus(tree_file), error = function(e) e)
  if (!inherits(tree_try, "error")) {
    return(if (inherits(tree_try, "multiPhylo")) tree_try[[1]] else tree_try)
  }

  message("ape::read.nexus failed; using fallback parser for MCC tree.")
  lines <- readLines(tree_file, warn = FALSE)

  tree_idx <- grep("^[[:space:]]*;?[Tt][Rr][Ee][Ee][[:space:]]+", lines)
  if (length(tree_idx) == 0) {
    stop("Could not find a tree line in MCC file: ", tree_file)
  }

  translate_idx <- grep("^[[:space:]]*[Tt]ranslate[[:space:]]*$", lines)
  translate_map <- NULL
  if (length(translate_idx) > 0 && translate_idx[1] < tree_idx[1]) {
    trans_lines <- trimws(lines[(translate_idx[1] + 1):(tree_idx[1] - 1)])
    trans_lines <- trans_lines[nzchar(trans_lines)]
    trans_lines <- gsub("[,;]$", "", trans_lines)

    parsed <- lapply(trans_lines, function(z) {
      m <- regexec("^([0-9]+)[[:space:]]+(.+)$", z)
      v <- regmatches(z, m)[[1]]
      if (length(v) != 3) return(NULL)
      c(id = v[2], label = trimws(v[3]))
    })
    parsed <- parsed[!vapply(parsed, is.null, logical(1))]
    if (length(parsed) > 0) {
      ids <- vapply(parsed, function(z) z["id"], character(1))
      labs <- vapply(parsed, function(z) z["label"], character(1))
      translate_map <- setNames(labs, ids)
    }
  }

  tree_start <- tree_idx[1]
  tail_lines <- lines[tree_start:length(lines)]
  end_rel <- grep(";", tail_lines)[1]
  if (is.na(end_rel)) {
    stop("Could not locate semicolon ending tree definition in: ", tree_file)
  }

  tree_lines <- tail_lines[seq_len(end_rel)]
  tree_text <- paste(tree_lines, collapse = " ")
  tree_text <- sub("^[[:space:]]*;?[Tt][Rr][Ee][Ee][[:space:]]+[^=]+=[[:space:]]*", "", tree_text)
  tree_text <- sub("[[:space:]]*[Ee][Nn][Dd];[[:space:]]*$", "", tree_text)
  tree_text <- trimws(tree_text)
  if (!grepl(";$", tree_text)) tree_text <- paste0(tree_text, ";")

  tree_obj <- ape::read.tree(text = tree_text)
  if (inherits(tree_obj, "multiPhylo")) {
    tip_counts <- vapply(tree_obj, ape::Ntip, integer(1))
    tree_obj <- tree_obj[[which.max(tip_counts)]]
  }

  if (!is.null(translate_map)) {
    mapped <- unname(translate_map[tree_obj$tip.label])
    keep_orig <- is.na(mapped)
    mapped[keep_orig] <- tree_obj$tip.label[keep_orig]
    tree_obj$tip.label <- mapped
  }

  tree_obj
}

compute_equal_splits <- function(tree) {
  es <- picante::evol.distinct(tree, type = "equal.splits")
  names(es)[1:2] <- c("glottocode", "ES")
  es %>%
    dplyr::transmute(
      glottocode = normalize_glottocode(glottocode),
      ES = as.numeric(ES)
    )
}

load_tip_metrics <- function(tree, ed_file) {
  fallback <- compute_equal_splits(tree) %>%
    dplyr::transmute(
      glottocode,
      ED_tip = ES,
      DR_tip = ifelse(ES > 0, 1 / ES, NA_real_)
    )

  ed_source <- "single-tree equal-splits from MCC (fallback)"
  dr_source <- "single-tree DR = 1/equal-splits from MCC (fallback)"
  metric_tbl <- NULL

  if (file.exists(ed_file)) {
    all_ed <- tryCatch(
      readr::read_csv(ed_file, show_col_types = FALSE, progress = FALSE),
      error = function(e) {
        warning("Could not read ALL_ED_wTREES.csv: ", e$message)
        NULL
      }
    )

    if (!is.null(all_ed)) {
      gc_col <- first_present_col(all_ed, c("glottocode", "Group.1", "label"))
      ed_col <- first_present_col(all_ed, c("ED", "w", "ED_tip"))

      if (!is.na(gc_col) && !is.na(ed_col)) {
        metric_tbl <- all_ed %>%
          dplyr::transmute(
            glottocode = normalize_glottocode(.data[[gc_col]]),
            ED = as.numeric(.data[[ed_col]])
          ) %>%
          dplyr::group_by(glottocode) %>%
          dplyr::summarise(ED_tip = mean(ED, na.rm = TRUE), .groups = "drop") %>%
          dplyr::mutate(DR_tip = ifelse(ED_tip > 0, 1 / ED_tip, NA_real_))

        ed_source <- "posterior mean ED from outputs/EDGEscores/ALL_ED_wTREES.csv"
        dr_source <- "posterior mean DR (1/ED) from outputs/EDGEscores/ALL_ED_wTREES.csv"
      } else {
        warning("ALL_ED_wTREES.csv found, but required columns were not detected. Using MCC fallback.")
      }
    }
  } else {
    warning("ALL_ED_wTREES.csv not found. Using MCC equal-splits fallback for ED and DR.")
  }

  if (is.null(metric_tbl)) {
    metric_tbl <- fallback
  }

  metric_tbl <- dplyr::tibble(glottocode = tree$tip.label) %>%
    dplyr::left_join(metric_tbl, by = "glottocode") %>%
    dplyr::left_join(fallback, by = "glottocode", suffix = c("", "_fallback")) %>%
    dplyr::mutate(
      ED_tip = ifelse(is.na(ED_tip), ED_tip_fallback, ED_tip),
      DR_tip = ifelse(is.na(DR_tip), DR_tip_fallback, DR_tip)
    ) %>%
    dplyr::select(glottocode, ED_tip, DR_tip)

  list(metrics = metric_tbl, ed_source = ed_source, dr_source = dr_source)
}

build_branch_metric <- function(tree, tip_metric_named, out_name) {
  n_tip <- ape::Ntip(tree)
  edge_child <- tree$edge[, 2]
  tip_metric_by_index <- as.numeric(tip_metric_named[tree$tip.label])

  out <- rep(NA_real_, length(edge_child))
  tip_mask <- edge_child <= n_tip
  out[tip_mask] <- tip_metric_by_index[edge_child[tip_mask]]

  internal_nodes <- unique(edge_child[!tip_mask])
  for (node_id in internal_nodes) {
    desc <- phytools::getDescendants(tree, node_id)
    desc_tips <- desc[desc <= n_tip]
    if (length(desc_tips) == 0) next
    out[edge_child == node_id] <- mean(tip_metric_by_index[desc_tips], na.rm = TRUE)
  }

  branch_df <- data.frame(node = edge_child, value = out)
  names(branch_df)[2] <- out_name
  branch_df
}

load_top27_map <- function(top27_file) {
  top27_raw <- readr::read_csv(top27_file, show_col_types = FALSE)
  fam_col <- first_present_col(top27_raw, c("Family", "family", "group"))
  tips_col <- first_present_col(top27_raw, c("tips", "tip", "glottocodes"))

  if (is.na(fam_col) || is.na(tips_col)) {
    stop("top27families.csv does not contain expected columns.")
  }

  top27 <- top27_raw %>%
    dplyr::transmute(Family = .data[[fam_col]], tips = .data[[tips_col]])

  family_map <- top27 %>%
    tidyr::separate_rows(tips, sep = ",") %>%
    dplyr::transmute(
      glottocode = normalize_glottocode(tips),
      top_family = Family
    ) %>%
    dplyr::filter(glottocode != "") %>%
    dplyr::distinct(glottocode, .keep_all = TRUE)

  list(map = family_map, levels = unique(top27$Family))
}

build_major_family_nodes <- function(tree, family_map, top_n = 10, min_tips = 30) {
  fam_counts <- family_map %>%
    dplyr::filter(glottocode %in% tree$tip.label) %>%
    dplyr::count(top_family, sort = TRUE) %>%
    dplyr::filter(n >= min_tips) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(fam_counts) == 0) return(data.frame())

  rows <- list()
  for (i in seq_len(nrow(fam_counts))) {
    fam <- fam_counts$top_family[i]
    fam_tips <- intersect(tree$tip.label, family_map$glottocode[family_map$top_family == fam])
    if (length(fam_tips) < 2) next

    mono <- tryCatch(ape::is.monophyletic(tree, fam_tips), error = function(e) FALSE)
    if (!isTRUE(mono)) next

    node_id <- tryCatch(ape::getMRCA(tree, fam_tips), error = function(e) NA_integer_)
    if (is.na(node_id)) next

    rows[[length(rows) + 1]] <- data.frame(
      top_family = fam,
      node = node_id,
      n_tips = length(fam_tips),
      label = pretty_family_name(fam),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) return(data.frame())
  dplyr::bind_rows(rows)
}

build_family_ring_labels <- function(tip_df, families_to_label) {
  if (length(families_to_label) == 0) return(data.frame())

  df <- tip_df %>%
    dplyr::mutate(
      idx = dplyr::row_number(),
      fam_chr = as.character(top_family)
    )
  n_tips <- nrow(df)
  out <- list()

  for (fam in families_to_label) {
    idxs <- df$idx[df$fam_chr == fam]
    if (length(idxs) == 0) next

    run_breaks <- c(1, which(diff(idxs) != 1) + 1)
    run_ends <- c(run_breaks[-1] - 1, length(idxs))
    run_lengths <- run_ends - run_breaks + 1
    best_run <- which.max(run_lengths)
    chosen_run <- idxs[run_breaks[best_run]:run_ends[best_run]]
    mid_idx <- chosen_run[ceiling(length(chosen_run) / 2)]

    theta <- 360 * ((mid_idx - 1) / n_tips)
    text_angle <- theta - 90
    text_hjust <- ifelse(text_angle < -90, 1, 0)
    if (text_angle < -90) text_angle <- text_angle + 180

    out[[length(out) + 1]] <- data.frame(
      label = df$label[mid_idx],
      ring_label = pretty_family_name(fam),
      text_angle = text_angle,
      text_hjust = text_hjust,
      stringsAsFactors = FALSE
    )
  }

  if (length(out) == 0) return(data.frame())
  dplyr::bind_rows(out)
}

load_predictor_table <- function(file_path, source_name, target_time3 = TARGET_TIME3) {
  if (!file.exists(file_path)) return(NULL)

  dat <- tryCatch(
    readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE),
    error = function(e) {
      warning("Could not read predictor file ", file_path, ": ", e$message)
      NULL
    }
  )
  if (is.null(dat)) return(NULL)

  # Prefer explicit glottocode columns over LANG_IS to avoid ID mismatches
  # in files that contain both.
  id_col <- first_present_col(dat, c("glottocode", "Glottocode", "LANG_IS"))
  if (is.na(id_col)) {
    warning("Skipping predictor file with no glottocode ID column: ", file_path)
    return(NULL)
  }

  if ("Time3" %in% names(dat)) {
    time_num <- suppressWarnings(as.numeric(dat$Time3))
    keep <- is.finite(time_num) & time_num == target_time3
    if (!any(keep)) {
      warning(
        "No rows with Time3 == ", target_time3,
        " in predictor file: ", file_path, ". Skipping this source."
      )
      return(NULL)
    }
    dat <- dat[keep, , drop = FALSE]
  }

  dplyr::tibble(
    glottocode = normalize_glottocode(dat[[id_col]]),
    area = pull_numeric_col(dat, c("Shap_Ar", "area", "Area", "Range_size", "range_size")),
    popd = pull_numeric_col(dat, c("popd", "population_density", "pop_density")),
    dist_city = pull_numeric_col(dat, c("distancetocityyear", "distancetocityyear2", "distance_to_city_year")),
    friction = pull_numeric_col(dat, c("friction", "mean_friction", "mean_friction_walking")),
    island = pull_numeric_col(dat, c("island", "Island")),
    source = source_name
  )
}

load_predictors <- function(tree_labels, path_in, path_out, target_time3 = TARGET_TIME3) {
  candidate_sources <- list(
    # Area-rich source used in environmental preprocessing (contains Shap_Ar)
    list(path = path_out("fixed", "glotto_languages_cites_states3.csv"), source = "outputs/fixed/glotto_languages_cites_states3.csv"),
    list(path = path_in("fixed", "final_env_datapoints.csv"), source = "input_data/fixed/final_env_datapoints.csv"),
    list(path = path_in("final_env_datapoints.csv"), source = "input_data/final_env_datapoints.csv"),
    list(path = path_out("fixed", "glotto_languages_cites_states3points.csv"), source = "outputs/fixed/glotto_languages_cites_states3points.csv"),
    list(path = path_in("bisse_data_in.csv"), source = "input_data/bisse_data_in.csv")
  )

  pred_list <- lapply(candidate_sources, function(x) {
    load_predictor_table(x$path, x$source, target_time3 = target_time3)
  })
  pred_list <- pred_list[!vapply(pred_list, is.null, logical(1))]

  if (length(pred_list) == 0) {
    empty <- dplyr::tibble(glottocode = tree_labels) %>%
      dplyr::mutate(
        area_raw = NA_real_,
        popd_raw = NA_real_,
        dist_city_raw = NA_real_,
        friction_raw = NA_real_,
        island_raw = NA_real_,
        area_log = NA_real_,
        popd_log = NA_real_,
        dist_city_log = NA_real_,
        friction_log = NA_real_,
        island_bin = factor(NA_character_, levels = c("Mainland", "Island"))
      )
    return(list(predictors = empty, sources = character(0)))
  }

  pred_raw <- dplyr::bind_rows(pred_list) %>%
    dplyr::filter(glottocode %in% tree_labels)

  pred_summary <- pred_raw %>%
    dplyr::group_by(glottocode) %>%
    dplyr::summarise(
      area = first_non_na(area),
      popd = first_non_na(popd),
      dist_city = first_non_na(dist_city),
      friction = first_non_na(friction),
      island = first_non_na(island),
      .groups = "drop"
    )

  pred <- dplyr::tibble(glottocode = tree_labels) %>%
    dplyr::left_join(pred_summary, by = "glottocode") %>%
    dplyr::mutate(
      area_raw = area,
      popd_raw = popd,
      dist_city_raw = dist_city,
      friction_raw = friction,
      island_raw = island,
      area_log = ifelse(is.na(area_raw), NA_real_, log1p(pmax(area_raw, 0))),
      popd_log = ifelse(is.na(popd_raw), NA_real_, log1p(pmax(popd_raw, 0))),
      dist_city_log = ifelse(is.na(dist_city_raw), NA_real_, log1p(pmax(dist_city_raw, 0))),
      friction_log = ifelse(is.na(friction_raw), NA_real_, log1p(pmax(friction_raw, 0))),
      island_bin = dplyr::case_when(
        is.na(island_raw) ~ NA_character_,
        island_raw > 0.5 ~ "Island",
        TRUE ~ "Mainland"
      ),
      island_bin = factor(island_bin, levels = c("Mainland", "Island"))
    ) %>%
    dplyr::select(
      glottocode,
      area_raw, popd_raw, dist_city_raw, friction_raw, island_raw,
      area_log, popd_log, dist_city_log, friction_log, island_bin
    )

  list(predictors = pred, sources = unique(pred_raw$source))
}

prepare_mcc_tree_variant_data <- function(
  base_dir = resolve_base_dir(),
  target_time3 = TARGET_TIME3,
  predictor_coverage_threshold = 0.80
) {
  ensure_mcc_data_packages()

  path_in <- function(...) file.path(base_dir, "input_data", ...)
  path_out <- function(...) file.path(base_dir, "outputs", ...)

  tree_file <- path_in("edge6636-March-2023-relabelled_cartoonised.tree")
  top27_file <- path_in("top27families.csv")
  ed_file <- path_out("EDGEscores", "ALL_ED_wTREES.csv")

  if (!file.exists(tree_file)) stop("MCC tree file not found: ", tree_file)
  if (!file.exists(top27_file)) stop("Top-27 family file not found: ", top27_file)

  tree <- read_mcc_tree_safe(tree_file)
  tree$tip.label <- sub("_.*$", "", tree$tip.label)
  tree$tip.label <- normalize_glottocode(tree$tip.label)
  if (anyDuplicated(tree$tip.label) > 0) {
    stop("Duplicate tip labels found after glottocode normalization.")
  }
  tree <- ape::ladderize(tree, right = FALSE)
  message("Loaded MCC tree with ", ape::Ntip(tree), " tips.")

  top27_info <- load_top27_map(top27_file)
  family_map <- top27_info$map
  family_levels <- top27_info$levels

  metric_info <- load_tip_metrics(tree, ed_file)
  tip_metrics <- metric_info$metrics

  tip_df <- dplyr::tibble(label = tree$tip.label) %>%
    dplyr::left_join(tip_metrics, by = c("label" = "glottocode")) %>%
    dplyr::left_join(family_map, by = c("label" = "glottocode")) %>%
    dplyr::mutate(
      top_family = ifelse(is.na(top_family), "Other", top_family),
      top_family = factor(top_family, levels = c(family_levels, "Other"))
    )

  if (nrow(tip_df) != ape::Ntip(tree)) {
    stop("Tip metadata row count does not match tree tip count.")
  }
  if (any(is.na(tip_df$ED_tip))) {
    stop("Missing ED values after fallback handling.")
  }

  ed_tip_clip <- clip_quantiles(tip_df$ED_tip, probs = c(0.01, 0.99))
  dr_tip_clip <- clip_quantiles(tip_df$DR_tip, probs = c(0.01, 0.99))
  tip_df$ED_tip_clip <- ed_tip_clip$values
  tip_df$DR_tip_clip <- dr_tip_clip$values

  tip_ed <- tip_df$ED_tip
  names(tip_ed) <- tip_df$label
  tip_dr <- tip_df$DR_tip
  names(tip_dr) <- tip_df$label

  branch_ed <- build_branch_metric(tree, tip_ed, "branch_ED")
  branch_dr <- build_branch_metric(tree, tip_dr, "branch_DR")

  ed_branch_clip <- clip_quantiles(branch_ed$branch_ED, probs = c(0.01, 0.99))
  dr_branch_clip <- clip_quantiles(branch_dr$branch_DR, probs = c(0.01, 0.99))
  branch_ed$branch_ED_clip <- ed_branch_clip$values
  branch_dr$branch_DR_clip <- dr_branch_clip$values

  family_nodes <- build_major_family_nodes(tree, family_map, top_n = 8, min_tips = 60)
  family_ring_labels <- build_family_ring_labels(tip_df, family_nodes$top_family)

  pred_info <- load_predictors(tree$tip.label, path_in, path_out, target_time3 = target_time3)
  tip_df <- tip_df %>%
    dplyr::left_join(pred_info$predictors, by = c("label" = "glottocode"))

  area_cov <- coverage_prop(tip_df$area_log)
  if (area_cov == 0) {
    stop(
      "Area data missing for all tips after predictor merge. ",
      "Provide one of: outputs/fixed/glotto_languages_cites_states3.csv ",
      "or input_data/bisse_data_in.csv."
    )
  }
  if (area_cov < predictor_coverage_threshold) {
    warning(
      "Area coverage below threshold (coverage = ", round(area_cov, 3),
      ", threshold = ", round(predictor_coverage_threshold, 3),
      "). Check glottocode matching across area sources."
    )
  }

  pred_coverage <- c(
    area_log = coverage_prop(tip_df$area_log),
    popd_log = coverage_prop(tip_df$popd_log),
    dist_city_log = coverage_prop(tip_df$dist_city_log),
    friction_log = coverage_prop(tip_df$friction_log)
  )

  eligible_cont <- names(pred_coverage[pred_coverage >= predictor_coverage_threshold])
  eligible_cont <- eligible_cont[order(pred_coverage[eligible_cont], decreasing = TRUE)]
  selected_cont <- head(eligible_cont, 4)

  show_island <- coverage_prop(tip_df$island_bin) >= predictor_coverage_threshold
  show_dr_ring <- length(selected_cont) == 0

  family_palette <- setNames(
    grDevices::hcl.colors(length(family_levels), palette = "Dark 3"),
    family_levels
  )
  family_palette <- c(family_palette, Other = "#E6E6E6")

  family_labels <- c(
    setNames(pretty_family_name(family_levels), family_levels),
    Other = "Other"
  )

  list(
    base_dir = base_dir,
    target_time3 = target_time3,
    predictor_coverage_threshold = predictor_coverage_threshold,
    paths = list(
      tree_file = tree_file,
      top27_file = top27_file,
      ed_file = ed_file,
      out_dir = path_out("Fig_MCC_tree_variants")
    ),
    tree = tree,
    tip_df = tip_df,
    branch_ed = branch_ed,
    branch_dr = branch_dr,
    ed_branch_clip = ed_branch_clip,
    dr_branch_clip = dr_branch_clip,
    ed_tip_clip = ed_tip_clip,
    dr_tip_clip = dr_tip_clip,
    family_nodes = family_nodes,
    family_ring_labels = family_ring_labels,
    family_palette = family_palette,
    family_labels = family_labels,
    pred_info = pred_info,
    pred_coverage = pred_coverage,
    selected_cont = selected_cont,
    show_island = show_island,
    show_dr_ring = show_dr_ring,
    metric_info = metric_info
  )
}
