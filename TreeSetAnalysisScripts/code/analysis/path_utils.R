ts_here <- function(...) {
  here::here("TreeSetAnalysisScripts", ...)
}

ts_source <- function(...) {
  source(ts_here(...), local = parent.frame())
}

ts_dir_create <- function(...) {
  dir_path <- ts_here(...)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  dir_path
}
