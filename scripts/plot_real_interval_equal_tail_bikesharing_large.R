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

read_metrics <- function(path_vec) {
  dfs <- lapply(path_vec, function(p) utils::read.csv(p, stringsAsFactors = FALSE))
  do.call(rbind, dfs)
}

prepare_plot_data <- function(df, dataset, type_levels, level_levels, allowed_models) {
  df <- df[df$Dataset == dataset & df$Model %in% allowed_models, , drop = FALSE]
  if (!nrow(df)) return(NULL)

  width_cols <- grep("^Width_(Q|HPD)_", names(df), value = TRUE)
  coverage_cols <- grep("^Coverage_(Q|HPD)_", names(df), value = TRUE)
  if (!length(width_cols) || !length(coverage_cols)) return(NULL)

  width_long <- tidyr::pivot_longer(
    df,
    cols = dplyr::all_of(width_cols),
    names_to = c("Type", "Level"),
    names_pattern = "Width_(Q|HPD)_(\\d+)",
    values_to = "Width",
    values_drop_na = TRUE
  )
  cov_long <- tidyr::pivot_longer(
    df,
    cols = dplyr::all_of(coverage_cols),
    names_to = c("Type", "Level"),
    names_pattern = "Coverage_(Q|HPD)_(\\d+)",
    values_to = "Coverage",
    values_drop_na = TRUE
  )

  keep_width <- width_long$Type %in% type_levels & width_long$Level %in% level_levels
  keep_cov <- cov_long$Type %in% type_levels & cov_long$Level %in% level_levels
  width_long <- width_long[keep_width, , drop = FALSE]
  cov_long <- cov_long[keep_cov, , drop = FALSE]
  if (!nrow(width_long) || !nrow(cov_long)) return(NULL)

  type_labels <- c(Q = "Equal-tail", HPD = "HPD")
  width_long$Type_Label <- type_labels[width_long$Type]
  width_long$Level_Label <- paste0(width_long$Level, "%")
  width_long$Model <- factor(width_long$Model, levels = allowed_models)

  cov_long$Type_Label <- type_labels[cov_long$Type]
  cov_long$Level_Label <- paste0(cov_long$Level, "%")
  cov_long$Model <- factor(cov_long$Model, levels = allowed_models)

  cov_ann <- dplyr::summarise(
    dplyr::group_by(cov_long, Model, Type_Label, Level_Label),
    Coverage = mean(Coverage, na.rm = TRUE),
    .groups = "drop"
  )
  width_max <- dplyr::summarise(
    dplyr::group_by(width_long, Model, Type_Label, Level_Label),
    MaxWidth = max(Width, na.rm = TRUE),
    .groups = "drop"
  )
  list(width = width_long, coverage = dplyr::left_join(cov_ann, width_max,
                                                        by = c("Model", "Type_Label", "Level_Label")))
}

build_plot <- function(width_df, cov_df, dataset_title, base_size = 16, label_size = 4.2) {
  fill_values <- c(
    sebart = "#1d4ed8",
    bart = "#9333ea",
    bart_bc = "#059669",
    sblm = "#f97316"
  )
  fill_values <- fill_values[names(fill_values) %in% levels(width_df$Model)]

  cov_df$Label <- sprintf("%s%%", round(cov_df$Coverage * 100))
  plot_obj <- ggplot2::ggplot(width_df, ggplot2::aes(x = Width, y = Model, fill = Model)) +
    ggplot2::geom_boxplot(alpha = 0.3, outlier.size = 0.9, linewidth = 0.4) +
    ggplot2::facet_grid(Type_Label ~ Level_Label, scales = "free_x") +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::geom_text(
      data = cov_df,
      ggplot2::aes(x = MaxWidth, y = Model, label = Label),
      colour = "#0f172a",
      hjust = -0.1,
      size = label_size
    ) +
    ggplot2::labs(
      title = sprintf("Interval widths: %s", dataset_title),
      x = "Interval width",
      y = NULL,
      fill = "Model"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey70"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text.x = ggplot2::element_text(size = base_size * 0.85, face = "bold"),
      strip.text.y = ggplot2::element_text(size = base_size * 0.85, face = "bold"),
      axis.text.y = ggplot2::element_text(size = base_size * 0.9, face = "bold"),
      axis.text.x = ggplot2::element_text(size = base_size * 0.8)
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.25)))
  plot_obj
}

main <- function() {
  root_dir <- resolve_root()
  args <- modifyList(list(
    metrics = file.path(root_dir, "results", "real_data", "real_forestfires_area_bikesharing_windspeed_parkinsons_DFA_metrics.csv"),
    dataset = "bikesharing",
    prefix = "real_interval_equal_tail",
    outdir = file.path(root_dir, "results", "real_data"),
    levels = "80,90,95",
    type = "Q",
    width = "9.5",
    height = "5.5",
    dpi = "450",
    base_size = "16"
  ), parse_args(commandArgs(trailingOnly = TRUE)))

  metric_paths <- strsplit(args$metrics, ",", fixed = TRUE)[[1]]
  metric_paths <- normalizePath(metric_paths, mustWork = TRUE)
  metrics_df <- read_metrics(metric_paths)

  dataset <- args$dataset
  levels_vec <- trimws(strsplit(args$levels, ",", fixed = TRUE)[[1]])
  type_vec <- trimws(strsplit(args$type, ",", fixed = TRUE)[[1]])
  allowed_models <- c("sebart", "bart", "bart_bc", "sblm")

  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
  })

  plot_data <- prepare_plot_data(metrics_df, dataset, type_vec, levels_vec, allowed_models)
  if (is.null(plot_data)) {
    stop(sprintf("No data available for dataset '%s' and the specified levels/types.", dataset))
  }

  dataset_title <- tools::toTitleCase(dataset)
  base_size <- as.numeric(args$base_size %||% 16)
  label_size <- base_size / 2.8
  plot_obj <- build_plot(plot_data$width, plot_data$coverage, dataset_title,
                         base_size = base_size, label_size = label_size)

  outdir <- normalizePath(args$outdir, mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(outdir, sprintf("%s_%s_interval_widths.png", args$prefix, dataset))
  ggplot2::ggsave(outfile, plot_obj,
                  width = as.numeric(args$width),
                  height = as.numeric(args$height),
                  dpi = as.numeric(args$dpi),
                  bg = "white")
  message("Saved plot to ", outfile)
}

if (sys.nframe() == 0) {
  main()
}
