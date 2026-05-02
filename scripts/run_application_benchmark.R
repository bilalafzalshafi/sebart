#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_application_benchmark.R --datasets=energy:heating_load,...

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

root_dir <- resolve_root()
options(sebart.root = root_dir)
source(file.path(root_dir, "sebart.R"))
source(file.path(root_dir, "bart_bc.R"))
source(file.path(root_dir, "R", "experiment_utils.R"))

run_application_benchmark_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  default_tasks <- paste(
    c(
      "energy:heating_load",
      "energy:cooling_load",
      "concrete:strength",
      "real_estate:price",
      "online_news:shares",
      "facebook_metrics:interactions",
      "power_consumption:zone_1",
      "power_consumption:zone_2",
      "power_consumption:zone_3",
      "stocks:total_risk",
      "stocks:systematic_risk",
      "stocks:annual_return",
      "productivity:actual_productivity",
      "math:math_grade",
      "portuguese:portuguese_grade",
      "bikesharing:windspeed",
      "forestfires:area",
      "parkinsons:DFA"
    ),
    collapse = ","
  )

  defaults <- list(
    datasets = default_tasks,
    models = "sebart,bart,bart_bc,sblm,drbart",
    folds = "10",
    repeats = "1",
    seed = "202405",
    alpha = "0.1",
    outdir = file.path(root_dir, "results", "real_data"),
    prefix = "application_benchmark",
    verbose = "TRUE"
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
  prefix <- opts$prefix
  verbose <- as_bool(opts$verbose)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  study <- run_real_data_study(
    dataset_targets = dataset_targets,
    models = models,
    folds = folds,
    repeats = repeats,
    seed = seed,
    alpha = alpha,
    verbose = verbose,
    ppc_models = intersect(models, c("sebart", "bart"))
  )

  metrics_path <- file.path(outdir, paste0(prefix, "_metrics.csv"))
  summary_path <- file.path(outdir, paste0(prefix, "_summary.csv"))
  write.csv(study$metrics, metrics_path, row.names = FALSE)
  write.csv(study$summary, summary_path, row.names = FALSE)

  relative_df <- compute_relative_benchmark_metrics(study$summary)
  if (!is.null(relative_df)) {
    rel_path <- file.path(outdir, paste0(prefix, "_relative_metrics.csv"))
    write.csv(relative_df, rel_path, row.names = FALSE)
    plots <- plot_relative_benchmark(relative_df)
    if (!is.null(plots) && length(plots)) {
      ensure_packages("ggplot2")
      plot_files <- list(
        `ΔRMSE` = file.path(outdir, paste0(prefix, "_relative_boxplot_delta_rmse.png")),
        `ΔCoverage Error` = file.path(outdir, paste0(prefix, "_relative_boxplot_delta_coverage.png")),
        `ΔWidth` = file.path(outdir, paste0(prefix, "_relative_boxplot_delta_width.png"))
      )
      for (plt in plots) {
        metric_name <- unique(plt$labels$title)
        file_path <- plot_files[[metric_name]]
        if (!is.null(file_path)) {
          ggplot2::ggsave(file_path, plt, width = 6, height = 4.5, dpi = 300, bg = "white")
          if (verbose) message("Saved benchmark plot to ", file_path)
        }
      }
    }
    if (verbose) message("Saved relative metrics to ", rel_path)
  }

  if (!is.null(study$metrics) && nrow(study$metrics)) {
    save_interval_widths_by_dataset(
      metrics_df = study$metrics,
      outdir = outdir,
      prefix = prefix,
      datasets = unique(study$metrics$Dataset)
    )
  }

  if (verbose) {
    message("Saved metrics to ", metrics_path)
    message("Saved summary to ", summary_path)
  }

  invisible(list(metrics_path = metrics_path, summary_path = summary_path))
}

if (sys.nframe() == 0) {
  run_application_benchmark_main()
}
