#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_simulations.R --iterations=100 --scenarios=sigmoid,beta

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

run_simulations_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    scenarios = paste(get_available_scenarios(), collapse = ","),
    iterations = "100",
    n_train = "200",
    n_test = "1000",
    p = "10",
    models = "sebart,bart,bart_bc,sblm",
    seed = "202401",
    predictor = "gaussian",
    alpha = "0.1",
    outdir = file.path(root_dir, "results", "simulations"),
    verbose = "TRUE",
    plot_iterations = "1"
  )

  parsed <- parse_args(args)
  opts <- modifyList(defaults, parsed)

  scenarios <- strsplit(opts$scenarios, ",")[[1]]
  iterations <- as.integer(opts$iterations)
  n_train <- as.integer(strsplit(opts$n_train, ",")[[1]])
  n_test <- as.integer(opts$n_test)
  p <- as.integer(opts$p)
  models <- strsplit(opts$models, ",")[[1]]
  seed <- as.integer(opts$seed)
  alpha <- as.numeric(opts$alpha)
  predictor <- opts$predictor
  outdir <- normalizePath(opts$outdir, mustWork = FALSE)
  verbose <- as_bool(opts$verbose)
  plot_iterations <- as.integer(opts$plot_iterations)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  study <- run_simulation_study(
    scenarios = scenarios,
    iterations = iterations,
    n_train = n_train,
    n_test = n_test,
    p = p,
    models = models,
    seed = seed,
    predictor = predictor,
    alpha = alpha,
    verbose = verbose,
    plot_iterations = plot_iterations
  )

  label <- paste0(
    "sim_",
    paste(scenarios, collapse = "-"),
    "_n",
    paste(n_train, collapse = "-")
  )
  label <- gsub("[^A-Za-z0-9]+", "_", label)

  metrics_path <- file.path(outdir, paste0(label, "_metrics.csv"))
  summary_path <- file.path(outdir, paste0(label, "_summary.csv"))
  write.csv(study$metrics, metrics_path, row.names = FALSE)
  write.csv(study$summary, summary_path, row.names = FALSE)

  if (!is.null(study$metrics) && nrow(study$metrics)) {
    n_values <- sort(unique(study$metrics$N_train))
    for (n_val in n_values) {
      prefix_n <- paste0(label, "_n", n_val)
      save_interval_widths_by_scenario(
        metrics_df = study$metrics,
        outdir = outdir,
        prefix = prefix_n,
        scenarios = unique(study$metrics$Scenario),
        n_train = n_val
      )
    }
  }

  if (!is.null(study$action_ranges)) {
    write.csv(study$action_ranges, file.path(outdir, paste0(label, "_latent_range.csv")), row.names = FALSE)
  }
  if (!is.null(study$transformation_metrics)) {
    write.csv(study$transformation_metrics, file.path(outdir, paste0(label, "_transformation_metrics.csv")), row.names = FALSE)
  }
  if (!is.null(study$transformation_curves)) {
    curve_csv <- file.path(outdir, paste0(label, "_transformation_curves.csv"))
    write.csv(study$transformation_curves, curve_csv, row.names = FALSE)
    save_transformation_plots(
      curve_df = study$transformation_curves,
      outdir = outdir,
      prefix = label,
      width = 10,
      height = 6,
      dpi = 300
    )
  }

  if (verbose) {
    message("Saved metrics to ", metrics_path)
    message("Saved summary to ", summary_path)
  }

  invisible(list(metrics_path = metrics_path, summary_path = summary_path))
}

if (sys.nframe() == 0) {
  run_simulations_main()
}
