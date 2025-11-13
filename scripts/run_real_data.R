#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_real_data.R --datasets=forestfires:area

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

root_dir <- resolve_root()
options(sebart.root = root_dir)
source(file.path(root_dir, "sebart.R"))
source(file.path(root_dir, "bart_bc.R"))
source(file.path(root_dir, "R", "experiment_utils.R"))
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

as_bool <- function(x) {
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

parse_dataset_targets <- function(spec) {
  entries <- strsplit(spec, ",")[[1]]
  lapply(entries, function(entry) {
    parts <- strsplit(entry, ":", fixed = TRUE)[[1]]
    dataset <- parts[1]
    available <- get_target_features(dataset)
    target <- if (length(parts) > 1) parts[2] else available[1]
    if (!target %in% available) {
      stop(sprintf("Target %s not available for dataset %s", target, dataset), call. = FALSE)
    }
    list(dataset = dataset, target = target)
  })
}

run_real_data_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    datasets = "forestfires:area,bikesharing:windspeed,parkinsons:DFA",
    models = "sebart,bart,bart_bc,sblm,lm",
    folds = "10",
    repeats = "10",
    seed = "202402",
    alpha = "0.1",
    outdir = file.path(root_dir, "results", "real_data"),
    verbose = "TRUE",
    ppc = "TRUE"
  )

  parsed <- parse_args(args)
  opts <- modifyList(defaults, parsed)

  dataset_targets <- parse_dataset_targets(opts$datasets)
  models <- strsplit(opts$models, ",")[[1]]
  folds <- as.integer(opts$folds)
  repeats <- as.integer(opts$repeats)
  seed <- as.integer(opts$seed)
  alpha <- as.numeric(opts$alpha)
  outdir <- normalizePath(opts$outdir, mustWork = FALSE)
  verbose <- as_bool(opts$verbose)
  make_ppc <- as_bool(opts$ppc)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  study <- run_real_data_study(
    dataset_targets = dataset_targets,
    models = models,
    folds = folds,
    repeats = repeats,
    seed = seed,
    alpha = alpha,
    verbose = verbose
  )

  label <- paste0("real_", paste(vapply(dataset_targets, function(dt) paste(dt$dataset, dt$target, sep = "-"), character(1)), collapse = "_"))
  label <- gsub("[^A-Za-z0-9]+", "_", label)

  metrics_path <- file.path(outdir, paste0(label, "_metrics.csv"))
  summary_path <- file.path(outdir, paste0(label, "_summary.csv"))
  write.csv(study$metrics, metrics_path, row.names = FALSE)
  write.csv(study$summary, summary_path, row.names = FALSE)

  dataset_tables <- create_dataset_tables(study$summary)
  lapply(names(dataset_tables), function(name) {
    table_path <- file.path(outdir, paste0(label, "_summary_", gsub("[^A-Za-z0-9]+", "_", name), ".csv"))
    write.csv(dataset_tables[[name]], table_path, row.names = FALSE)
  })

  width_plot <- plot_real_data_widths(study$metrics, title = "Real-data interval widths")
  if (!is.null(width_plot)) {
    plot_path <- file.path(outdir, paste0(label, "_interval_widths.png"))
    ensure_packages("ggplot2")
    ggplot2::ggsave(plot_path, width_plot, width = 8, height = 5, dpi = 300, bg = "white")
  }

  if (!is.null(study$metrics) && nrow(study$metrics)) {
    save_interval_widths_by_dataset(
      metrics_df = study$metrics,
      outdir = outdir,
      prefix = label,
      datasets = unique(study$metrics$Dataset)
    )
  }

  if (!is.null(study$posterior_samples)) {
    write.csv(study$posterior_samples, file.path(outdir, paste0(label, "_posterior_samples.csv")), row.names = FALSE)
  }
  if (!is.null(study$observed)) {
    write.csv(study$observed, file.path(outdir, paste0(label, "_observed.csv")), row.names = FALSE)
  }

  if (make_ppc && !is.null(study$posterior_samples) && !is.null(study$observed)) {
    ensure_packages("ggplot2")
    combos <- unique(study$posterior_samples[, c("Dataset", "Target")])
    apply(combos, 1, function(row) {
      dataset <- row[["Dataset"]]
      target <- row[["Target"]]
      ecdf_plot <- plot_ecdf_ppc(study$posterior_samples, study$observed, dataset, target)
      if (!is.null(ecdf_plot)) {
        out_file <- file.path(outdir, paste0(label, "_posterior_ecdf_", dataset, "_", target, ".png"))
        ggplot2::ggsave(out_file, ecdf_plot, width = 8, height = 6, dpi = 300, bg = "white")
      }
      density_plot <- plot_density_ppc(study$posterior_samples, study$observed, dataset, target)
      if (!is.null(density_plot)) {
        out_file <- file.path(outdir, paste0(label, "_posterior_density_", dataset, "_", target, ".png"))
        ggplot2::ggsave(out_file, density_plot, width = 8, height = 5, dpi = 300, bg = "white")
      }
      NULL
    })
  }

  if (verbose) {
    message("Saved metrics to ", metrics_path)
    message("Saved summary to ", summary_path)
  }

  invisible(list(metrics_path = metrics_path, summary_path = summary_path))
}

if (sys.nframe() == 0) {
  run_real_data_main()
}
