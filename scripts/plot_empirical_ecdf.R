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
  root <- resolve_root()
  options(sebart.root = root)
  source(file.path(root, "R", "experiment_utils.R"))

  defaults <- list(
    dataset = "parkinsons",
    target = "DFA",
    outdir = file.path(root, "results", "real_data"),
    prefix = "empirical_ecdf"
  )

  args <- modifyList(defaults, parse_args(commandArgs(trailingOnly = TRUE)))
  dataset_name <- args$dataset
  target_feature <- args$target
  outdir <- normalizePath(args$outdir, mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prefix <- args$prefix %||% "empirical_ecdf"

  data_obj <- load_real_dataset(dataset_name, target_feature, scale_y = FALSE, scale_X = FALSE)
  y_vec <- data_obj$y
  if (!length(y_vec)) stop("No observations found for dataset/target combination", call. = FALSE)

  ensure_packages("ggplot2")
  plot_df <- data.frame(Value = y_vec)
  title_txt <- sprintf("Empirical CDF: %s (%s)", tools::toTitleCase(dataset_name), target_feature)
  plt_ecdf <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Value)) +
    ggplot2::stat_ecdf(geom = "step", colour = "black", linewidth = 1) +
    ggplot2::labs(title = title_txt, x = "Value", y = "ECDF") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  outfile_ecdf <- file.path(outdir, sprintf("%s_%s_%s.png", prefix, dataset_name, target_feature))
  ggplot2::ggsave(outfile_ecdf, plt_ecdf, width = 7.5, height = 5.5, dpi = 300, bg = "white")
  message("Saved ECDF plot to ", outfile_ecdf)

  dens <- stats::density(y_vec)
  dens_df <- data.frame(x = dens$x, y = dens$y)
  title_density <- sprintf("Empirical Density: %s (%s)", tools::toTitleCase(dataset_name), target_feature)
  plt_density <- ggplot2::ggplot(dens_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(colour = "#2563eb", linewidth = 1.1) +
    ggplot2::labs(title = title_density, x = "Value", y = "Density") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  outfile_density <- file.path(outdir, sprintf("%s_%s_%s_density.png", prefix, dataset_name, target_feature))
  ggplot2::ggsave(outfile_density, plt_density, width = 7.5, height = 5.5, dpi = 300, bg = "white")
  message("Saved density plot to ", outfile_density)
}

if (sys.nframe() == 0) {
  main()
}
