locate_project_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  current <- start
  repeat {
    markers <- c("app.R", "helpers.R", "data")
    if (all(file.exists(file.path(current, markers)))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate project root from: ", start)
    }
    current <- parent
  }
}

init_project <- function(knitr_root = FALSE) {
  root <- locate_project_root()
  if (knitr_root && requireNamespace("knitr", quietly = TRUE)) {
    knitr::opts_knit$set(root.dir = root)
  }

  paths <- list(
    root = root,
    data_dir = file.path(root, "data"),
    plots_dir = file.path(root, "plots"),
    results_dir = file.path(root, "results"),
    pca_dir = file.path(root, "results", "pca"),
    wgcna_dir = file.path(root, "results", "wgcna"),
    functional_dir = file.path(root, "results", "functional_analysis"),
    gsea_dir = file.path(root, "results", "functional_analysis_gsea"),
    ora_dir = file.path(root, "results", "functional_analysis_ora"),
    network_dir = file.path(root, "results", "network_annotation"),
    supplementary_dir = file.path(root, "results", "supplementary")
  )

  dirs_to_create <- paths[c(
    "plots_dir",
    "results_dir",
    "pca_dir",
    "wgcna_dir",
    "functional_dir",
    "gsea_dir",
    "ora_dir",
    "network_dir",
    "supplementary_dir"
  )]

  invisible(lapply(dirs_to_create, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}
