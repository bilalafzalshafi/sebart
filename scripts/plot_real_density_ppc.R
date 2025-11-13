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
    samples = file.path(root_dir, "results", "real_data", "real_forestfires_area_bikesharing_windspeed_parkinsons_DFA_posterior_samples.csv"),
    datasets = "forestfires:area,bikesharing:windspeed,parkinsons:DFA",
    outdir = file.path(root_dir, "results", "real_data"),
    prefix = "real_density"
  )

  args <- modifyList(defaults, parse_args(commandArgs(trailingOnly = TRUE)))
  sample_paths <- strsplit(args$samples, ",")[[1]]
  sample_paths <- normalizePath(sample_paths, mustWork = TRUE)
  samples_list <- lapply(sample_paths, function(p) utils::read.csv(p, stringsAsFactors = FALSE))
  samples_df <- do.call(rbind, samples_list)
  samples_df <- samples_df[!is.na(samples_df$Model) & nzchar(samples_df$Model), , drop = FALSE]

  dataset_specs <- strsplit(args$datasets, ",")[[1]]
  dataset_specs <- trimws(dataset_specs)
  dataset_targets <- do.call(rbind, lapply(dataset_specs, function(spec) {
    parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
    data.frame(Dataset = parts[1], Target = parts[if (length(parts) > 1) 2 else 1], stringsAsFactors = FALSE)
  }))

  observed_list <- lapply(seq_len(nrow(dataset_targets)), function(i) {
    ds <- dataset_targets$Dataset[i]
    target <- dataset_targets$Target[i]
    dat <- load_real_dataset(ds, target)
    data.frame(
      Dataset = ds,
      Target = target,
      Observed = dat$y,
      stringsAsFactors = FALSE
    )
  })
  observed_df <- do.call(rbind, observed_list)

  dataset_vec <- unique(dataset_targets$Dataset)
  outdir <- normalizePath(args$outdir, mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prefix <- args$prefix %||% "real_density"

  saved <- character()
  for (i in seq_len(nrow(dataset_targets))) {
    dataset <- dataset_targets$Dataset[i]
    target <- dataset_targets$Target[i]
    plt <- plot_density_ppc(samples_df, observed_df, dataset, target)
    if (is.null(plt)) next
    outfile <- file.path(outdir, sprintf("%s_%s_%s_density.png", prefix, dataset, target))
    ggplot2::ggsave(outfile, plt, width = 7, height = 4.5, dpi = 300, bg = "white")
    saved <- c(saved, outfile)
  }

  if (length(saved)) {
    message("Saved PPC density plots:\n", paste(saved, collapse = "\n"))
  } else {
    message("No density plots generated (check dataset/target filters).")
  }
}

if (sys.nframe() == 0) {
  main()
}
