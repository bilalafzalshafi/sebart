#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_real_data_interpretability.R --dataset=forestfires --target=area
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
compute_varcount_importance <- function(varcount, feature_names) {
  totals <- colSums(varcount, na.rm = TRUE)
  if (!length(totals)) return(NULL)
  data.frame(
    Feature = feature_names,
    Importance = totals / sum(totals),
    stringsAsFactors = FALSE
  )
}
predict_sebart_grid <- function(fit, X_grid) {
  if (is.null(fit$post_ypred)) return(NULL)
  matrix_stats <- apply(fit$post_ypred, 2, mean)
  data.frame(Value = matrix_stats)
}
# effect_summary_from_draws no longer used in this script; the real-data pipeline now handles effect testing.
run_interpretability_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  root_dir <- resolve_root()
  options(sebart.root = root_dir)
  source(file.path(root_dir, "R", "experiment_utils.R"))
  source(file.path(root_dir, "sebart.R"))
  source(file.path(root_dir, "bart_bc.R"))
  defaults <- list(
    dataset = "forestfires",
    target = NA_character_,
    seed = "2025",
    outdir = file.path(root_dir, "results", "interpretability"),
    top_variables = "3",
    grid_points = "50",
    verbose = "TRUE"
  )
  opts <- modifyList(defaults, parse_args(args))
  dataset_name <- opts$dataset
  target_feature <- if (!is.na(opts$target)) opts$target else NULL
  seed <- as.integer(opts$seed)
  outdir <- normalizePath(opts$outdir, mustWork = FALSE)
  top_k <- max(1L, as.integer(opts$top_variables))
  grid_points <- max(10L, as.integer(opts$grid_points))
  verbose <- as_bool(opts$verbose)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ensure_packages(c("ggplot2"))
  if (verbose) message(sprintf("Loading dataset %s (%s)", dataset_name, target_feature %||% "default"))
  data <- load_real_dataset(dataset_name, target_feature)
  X <- data$X
  y <- data$y
  feature_names <- paste0("X", seq_len(ncol(X)))
  set.seed(seed)
  if (verbose) message("Fitting SeBART with variable selection...")
  sebart_fit <- sebart(
    y = y,
    X = X,
    X_test = X,
    ntree = 200,
    nsave = 1000,
    nburn = 1000,
    verbose = FALSE,
    var_select = TRUE,
    var_select_threshold = 0
  )
  if (verbose) message("Fitting standard BART...")
  bart_fit <- dbarts::bart(
    x.train = X,
    y.train = y,
    x.test = X,
    ntree = 200,
    ndpost = 1000,
    nskip = 1000,
    verbose = FALSE
  )
  var_importance <- list()
  if (!is.null(sebart_fit$var_importance_scores)) {
    var_importance$sebart <- data.frame(
      Feature = feature_names,
      Importance = sebart_fit$var_importance_scores,
      Model = "sebart",
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(bart_fit$varcount)) {
    lm_importance <- compute_varcount_importance(bart_fit$varcount, feature_names)
    if (!is.null(lm_importance)) {
      lm_importance$Model <- "bart"
      var_importance$bart <- lm_importance
    }
  }
  if (!is.null(var_importance) && length(var_importance)) {
    importance_df <- do.call(rbind, var_importance)
    write.csv(importance_df, file.path(outdir, sprintf("interpretability_%s_%s_importance.csv", dataset_name, target_feature %||% "default")), row.names = FALSE)
  }
  if (verbose) message("Computing partial dependence curves...")
  ranked_features <- NULL
  if (!is.null(sebart_fit$var_importance_scores)) {
    ranked_features <- order(sebart_fit$var_importance_scores, decreasing = TRUE)
  } else if (!is.null(var_importance$bart)) {
    ranked_features <- order(var_importance$bart$Importance, decreasing = TRUE)
  } else {
    ranked_features <- seq_len(min(3, ncol(X)))
  }
  ranked_features <- ranked_features[seq_len(min(top_k, length(ranked_features)))]
  pd_results <- list()
  for (idx in ranked_features) {
    grid_vals <- seq(min(X[, idx]), max(X[, idx]), length.out = grid_points)
    X_grid <- matrix(rep(colMeans(X), each = grid_points), nrow = grid_points)
    X_grid[, idx] <- grid_vals
    se_grid_fit <- sebart(
      y = y,
      X = X,
      X_test = X_grid,
      ntree = 200,
      nsave = 1000,
      nburn = 1000,
      verbose = FALSE
    )
    se_pd <- colMeans(se_grid_fit$post_ypred)
    bart_pred <- dbarts::predict(bart_fit, newdata = X_grid)
    bart_pd <- colMeans(bart_pred)
    pd_results[[length(pd_results) + 1]] <- data.frame(
      Feature = feature_names[idx],
      Value = grid_vals,
      Prediction = se_pd,
      Model = "sebart",
      stringsAsFactors = FALSE
    )
    pd_results[[length(pd_results) + 1]] <- data.frame(
      Feature = feature_names[idx],
      Value = grid_vals,
      Prediction = bart_pd,
      Model = "bart",
      stringsAsFactors = FALSE
    )
  }
  if (length(pd_results)) {
    pd_df <- do.call(rbind, pd_results)
    out_path <- file.path(outdir, sprintf("interpretability_%s_%s_partial_dependence.csv", dataset_name, target_feature %||% "default"))
    write.csv(pd_df, out_path, row.names = FALSE)
    plt <- ggplot2::ggplot(pd_df, ggplot2::aes(x = Value, y = Prediction, colour = Model)) +
    ggplot2::geom_line(size = 0.9) +
      ggplot2::facet_wrap(~Feature, scales = "free_x") +
      ggplot2::labs(title = sprintf("Partial dependence: %s - %s", dataset_name, target_feature %||% "default"), x = "Feature value", y = "Prediction") +
      ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(file.path(outdir, sprintf("interpretability_%s_%s_partial_dependence.png", dataset_name, target_feature %||% "default")), plt, width = 8, height = 5, dpi = 300, bg = "white")
  }

  if (verbose) message("Interpretability artefacts saved to ", outdir)
}
