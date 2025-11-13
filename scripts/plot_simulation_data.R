#!/usr/bin/env Rscript
# Usage: Rscript scripts/plot_simulation_data.R --scenarios=sigmoid --n_train=200

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

run_overlay_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  root_dir <- resolve_root()
  options(sebart.root = root_dir)
  source(file.path(root_dir, "R", "experiment_utils.R"))

  defaults <- list(
    scenarios = paste(get_available_scenarios(), collapse = ","),
    n_train = "200",
    n_test = "1000",
    p = "10",
    seed = "202401",
    predictor = "gaussian",
    outdir = file.path(root_dir, "results", "simulations"),
    verbose = "TRUE"
  )

  parsed <- parse_args(args)
  opts <- modifyList(defaults, parsed)

  scenarios <- strsplit(opts$scenarios, ",")[[1]]
  n_train <- as.integer(opts$n_train)
  n_test <- as.integer(opts$n_test)
  p <- as.integer(opts$p)
  seed <- as.integer(opts$seed)
  predictor <- opts$predictor
  outdir <- normalizePath(opts$outdir, mustWork = FALSE)
  verbose <- as_bool(opts$verbose)

  ensure_packages(c("ggplot2", "dplyr", "tidyr"))
  `%>%` <- dplyr::`%>%`
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  overlay_records <- list()

  for (idx in seq_along(scenarios)) {
    scenario <- scenarios[idx]
    if (verbose) message(sprintf("Generating overlay for scenario %s", scenario))
    data <- simulate_sebart_data(
      n_train = n_train,
      n_test = n_test,
      p = p,
      scenario = scenario,
      seed = seed + idx - 1,
      predictor = predictor
    )

    df <- dplyr::bind_rows(
      data.frame(value = data$y_train, split = "train", quantity = "Observed y"),
      data.frame(value = data$y_test, split = "test", quantity = "Observed y"),
      data.frame(value = data$z_train, split = "train", quantity = "Latent z"),
      data.frame(value = data$z_test, split = "test", quantity = "Latent z"),
      data.frame(value = data$f_true_train, split = "train", quantity = "f(x)"),
      data.frame(value = data$f_true_test, split = "test", quantity = "f(x)")
    )

    plot_title <- sprintf("Scenario %s", scenario)
    plt <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = split, colour = split)) +
      ggplot2::geom_density(alpha = 0.25) +
      ggplot2::facet_wrap(~quantity, scales = "free", ncol = 1) +
      ggplot2::labs(title = plot_title, x = NULL, y = "Density", fill = NULL, colour = NULL) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")

    out_file <- file.path(outdir, sprintf("sim_%s_overlay.png", gsub("[^A-Za-z0-9]+", "_", scenario)))
    ggplot2::ggsave(out_file, plt, width = 7, height = 8, dpi = 300, bg = "white")

    overlay_records[[scenario]] <- data.frame(
      Scenario = scenario,
      Split = rep(c("train", "test"), c(length(data$y_train), length(data$y_test))),
      Value = c(data$y_train, data$y_test),
      stringsAsFactors = FALSE
    )
  }

  if (length(overlay_records) > 1) {
    combined <- dplyr::bind_rows(overlay_records)
    combined$Scenario <- factor(combined$Scenario, levels = scenarios)

    combined <- combined %>%
      dplyr::group_by(Scenario, Split) %>%
      dplyr::mutate(
        Mean = mean(Value),
        Sd = stats::sd(Value)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Normalized = (Value - Mean) / (ifelse(Sd < 1e-8, 1e-8, Sd)))

    combined_plot <- ggplot2::ggplot(combined, ggplot2::aes(x = Normalized, colour = Scenario)) +
      ggplot2::geom_density(linewidth = 0.9, adjust = 1) +
      ggplot2::facet_wrap(~Split, ncol = 1, scales = "free") +
      ggplot2::labs(
        title = "Observed response overlays across scenarios (standardised)",
        x = "Standardised observed y",
        y = "Density",
        colour = "Scenario"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")

    combined_file <- file.path(outdir, paste0("sim_", paste(scenarios, collapse = "-"), "_observed_overlay.png"))
    ggplot2::ggsave(combined_file, combined_plot, width = 7, height = 6, dpi = 300, bg = "white")
  }
}

if (sys.nframe() == 0) {
  run_overlay_main()
}
