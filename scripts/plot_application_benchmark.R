#!/usr/bin/env Rscript
# Usage: Rscript scripts/plot_application_benchmark.R --summary=path --metrics=path

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

as_bool <- function(x) tolower(x) %in% c("true", "t", "1", "yes", "y")

parse_file_list <- function(value) {
  entries <- unlist(strsplit(value, ","))
  entries <- trimws(entries)
  entries <- entries[nzchar(entries)]
  vapply(entries, normalizePath, character(1), mustWork = TRUE)
}

flag_outliers <- function(values, multiplier = 6) {
  if (length(values) == 0) return(logical(0))
  center <- stats::median(values, na.rm = TRUE)
  scale <- stats::mad(values, center = center, constant = 1, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) {
    q <- stats::quantile(values, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    scale <- (q[2] - q[1]) / 1.349
  }
  if (!is.finite(scale) || scale == 0) {
    return(rep(FALSE, length(values)))
  }
  abs(values - center) > multiplier * scale
}

remove_metric_outliers <- function(df, metric_col, metric_label,
                                  multiplier = 6, max_remove = Inf) {
  if (!metric_col %in% names(df)) return(list(filtered = df, removed = df[FALSE, ]))
  values <- df[[metric_col]]
  mask <- flag_outliers(values, multiplier)
  removed <- df[mask, , drop = FALSE]
  if (!nrow(removed)) return(list(filtered = df, removed = removed))
  idx <- which(mask)
  order_idx <- order(abs(values[idx]), decreasing = TRUE)
  keep <- idx[order_idx[seq_len(min(length(order_idx), max_remove))]]
  removed <- df[keep, , drop = FALSE]
  removed$Metric <- metric_label
  removed$Value <- removed[[metric_col]]
  list(filtered = df[-keep, , drop = FALSE], removed = removed)
}

root_dir <- resolve_root()
options(sebart.root = root_dir)
source(file.path(root_dir, "R", "experiment_utils.R"))

plot_application_benchmark_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    summary = file.path(root_dir, "results", "real_data", "application_benchmark_summary.csv"),
    metrics = file.path(root_dir, "results", "real_data", "application_benchmark_metrics.csv"),
    outdir = file.path(root_dir, "results", "real_data"),
    prefix = "application_benchmark",
    verbose = "TRUE",
    filter_outliers = "TRUE",
    outlier_multiplier = "8",
    max_outliers_per_metric = "2"
  )

  opts <- modifyList(defaults, parse_args(args))
  verbose <- as_bool(opts$verbose)

  summary_paths <- parse_file_list(opts$summary)
  metrics_paths <- parse_file_list(opts$metrics)
  outdir <- normalizePath(opts$outdir, mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prefix <- opts$prefix
  filter_outliers <- as_bool(opts$filter_outliers)
  outlier_mult <- as.numeric(opts$outlier_multiplier)
  if (!is.finite(outlier_mult) || outlier_mult <= 0) outlier_mult <- 6
  max_outliers <- as.integer(opts$max_outliers_per_metric)
  if (!is.finite(max_outliers) || max_outliers <= 0) max_outliers <- 3

  summary_df <- do.call(rbind, lapply(summary_paths, utils::read.csv, stringsAsFactors = FALSE))
  if ("Repeat" %in% names(summary_df)) {
    summary_df <- summary_df[summary_df$Repeat == 1, , drop = FALSE]
  }

  metrics_df <- do.call(rbind, lapply(metrics_paths, utils::read.csv, stringsAsFactors = FALSE))
  if ("Repeat" %in% names(metrics_df)) {
    metrics_df <- metrics_df[metrics_df$Repeat == 1, , drop = FALSE]
  }

  relative_df <- compute_relative_benchmark_metrics(summary_df)
  if (!is.null(relative_df)) {
    rel_path <- file.path(outdir, paste0(prefix, "_relative_metrics.csv"))
    write.csv(relative_df, rel_path, row.names = FALSE)
    filtered_relative <- relative_df
    removed_df <- relative_df[FALSE, ]
    if (filter_outliers) {
      res_rmse <- remove_metric_outliers(
        filtered_relative, "Delta_RMSE", "ΔRMSE",
        multiplier = outlier_mult, max_remove = max_outliers
      )
      filtered_relative <- res_rmse$filtered
      removed_df <- rbind(removed_df, res_rmse$removed[, c("Dataset", "Target", "Model", "Metric", "Value")])

      res_width <- remove_metric_outliers(
        filtered_relative, "Delta_Width", "ΔWidth",
        multiplier = outlier_mult, max_remove = max_outliers
      )
      filtered_relative <- res_width$filtered
      removed_df <- rbind(removed_df, res_width$removed[, c("Dataset", "Target", "Model", "Metric", "Value")])

      if (nrow(removed_df)) {
        removed_path <- file.path(outdir, paste0(prefix, "_relative_outliers.csv"))
        write.csv(removed_df, removed_path, row.names = FALSE)
        if (verbose) {
          message("Removed ", nrow(removed_df), " outlier(s) before plotting. See ", removed_path)
        }
      }
    }
    plots <- plot_relative_benchmark(filtered_relative)
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
          ggplot2::ggsave(file_path, plt + ggplot2::theme(legend.position = "none"),
                          width = 6, height = 4.5, dpi = 300, bg = "white")
          if (verbose) message("Saved benchmark plot to ", file_path)
        }
      }
    }
    if (verbose) message("Saved relative metrics to ", rel_path)
  } else if (verbose) {
    message("No relative metrics generated (check summary file).")
  }

  if (!is.null(metrics_df) && nrow(metrics_df)) {
    datasets <- unique(metrics_df$Dataset)
    width_files <- save_interval_widths_by_dataset(
      metrics_df = metrics_df,
      outdir = outdir,
      prefix = prefix,
      datasets = datasets,
      levels = c("90"),
      interval_types = c("Q"),
      allowed_models = c("sebart", "drbart", "bart", "bart_bc", "sblm")
    )
    if (verbose && length(width_files)) {
      message("Saved interval-width plots:\n", paste(width_files, collapse = "\n"))
    }
  } else if (verbose) {
    message("Metrics file had no rows; skipping interval-width plots.")
  }

  invisible(list(summary = summary_paths, metrics = metrics_paths))
}

if (sys.nframe() == 0) {
  plot_application_benchmark_main()
}
