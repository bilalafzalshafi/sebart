#!/usr/bin/env Rscript

resolve_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) && nzchar(file_arg)) {
    script_dir <- dirname(normalizePath(file_arg))
  } else {
    ofile <- NULL
    frame1 <- sys.frames()[[1]]
    if (!is.null(frame1$ofile)) {
      ofile <- frame1$ofile
    } else if (!is.null(frame1$`__file__`)) {
      ofile <- frame1$`__file__`
    }
    script_dir <- if (!is.null(ofile)) dirname(normalizePath(ofile)) else normalizePath(getwd())
  }
  if (basename(script_dir) %in% c("scripts", "bin")) {
    normalizePath(file.path(script_dir, ".."))
  } else {
    script_dir
  }
}

parse_args <- function(args) {
  parsed <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    value <- if (length(kv) > 1) kv[2] else NA_character_
    parsed[[key]] <- value
  }
  parsed
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)) a else b

main <- function() {
  root_dir <- resolve_root()
  options(sebart.root = root_dir)
  source(file.path(root_dir, "R", "experiment_utils.R"))

  defaults <- list(
    metrics = file.path(root_dir, "results", "real_data", "real_forestfires_area_bikesharing_windspeed_parkinsons_DFA_metrics.csv"),
    datasets = "forestfires,bikesharing,parkinsons",
    outdir = file.path(root_dir, "results", "real_data"),
    prefix = "real_interval"
  )

  args <- modifyList(defaults, parse_args(commandArgs(trailingOnly = TRUE)))
  metric_paths <- strsplit(args$metrics, ",")[[1]]
  metric_paths <- normalizePath(metric_paths, mustWork = TRUE)
  metrics_list <- lapply(metric_paths, function(p) utils::read.csv(p, stringsAsFactors = FALSE))
  metrics_df <- do.call(rbind, metrics_list)

  dataset_vec <- strsplit(args$datasets, ",")[[1]]
  dataset_vec <- trimws(dataset_vec)
  outdir <- normalizePath(args$outdir, mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prefix <- args$prefix %||% "real_interval"

  saved <- save_interval_widths_by_dataset(
    metrics_df = metrics_df,
    outdir = outdir,
    prefix = prefix,
    datasets = dataset_vec
  )

  if (length(saved)) {
    message("Saved interval-width plots:\n", paste(saved, collapse = "\n"))
  } else {
    message("No plots generated (check dataset filter).")
  }
}

if (sys.nframe() == 0) {
  main()
}
