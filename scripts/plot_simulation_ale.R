#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dbarts)
})

friedman_function <- function(X) {
  X <- as.matrix(X)
  10 * sin(pi * X[, 1] * X[, 2]) +
    20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] +
    5 * X[, 5]
}

oracle_observed_mean <- function(X, transform_info, n_mc = 10000, seed = 123,
                                 chunk_size = 500) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  eps <- rnorm(n_mc)
  eps <- eps - mean(eps)

  f <- friedman_function(X)
  out <- numeric(length(f))

  for (start in seq(1, length(f), by = chunk_size)) {
    idx <- start:min(start + chunk_size - 1, length(f))
    z_mat <- outer(f[idx], eps, "+")
    y_mat <- transform_info$g_inv(z_mat)
    if (is.null(dim(y_mat))) {
      y_mat <- matrix(y_mat, nrow = length(idx), ncol = length(eps))
    }
    out[idx] <- rowMeans(y_mat)
  }

  out
}

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
as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")

jitter_bins <- function(edges) {
  edges <- as.numeric(edges)
  for (idx in seq_along(edges)) {
    if (idx > 1 && edges[idx] <= edges[idx - 1]) {
      edges[idx] <- edges[idx - 1] + 1e-6
    }
  }
  edges
}

build_ale_design <- function(X, feature_idx, K) {
  n <- nrow(X)
  block_list <- list(base = X)
  feature_info <- list()
  for (j in feature_idx) {
    xj <- X[, j]
    if (length(unique(xj)) < 2) next
    bins <- stats::quantile(xj, probs = seq(0, 1, length.out = K + 1), type = 1, na.rm = TRUE)
    bins[1] <- min(xj)
    bins[length(bins)] <- max(xj)
    bins <- jitter_bins(bins)
    block_names <- character(length(bins))
    for (idx in seq_along(bins)) {
      block_name <- sprintf("feat%s_b%s", j, idx)
      X_mut <- X
      X_mut[, j] <- bins[idx]
      block_list[[block_name]] <- X_mut
      block_names[idx] <- block_name
    }
    feature_info[[as.character(j)]] <- list(bins = bins, blocks = block_names)
  }
  block_sizes <- vapply(block_list, nrow, integer(1))
  cum_sizes <- cumsum(block_sizes)
  start_idx <- c(0, cum_sizes[-length(cum_sizes)])
  block_rows <- lapply(seq_along(block_list), function(idx) {
    if (block_sizes[idx] == 0) integer(0) else (start_idx[idx] + 1):cum_sizes[idx]
  })
  names(block_rows) <- names(block_list)

  X_test <- do.call(rbind, block_list)
  list(X_test = X_test,
       block_rows = block_rows,
       feature_info = feature_info,
       n = n)
}

compute_ale <- function(xj, bins, boundary_preds) {
  K <- length(bins) - 1
  if (K <= 0) return(NULL)
  ale_vals <- numeric(K + 1)
  for (k in seq_len(K)) {
    idx <- which(xj >= bins[k] & xj <= bins[k + 1])
    if (!length(idx)) {
      avg <- 0
    } else {
      diffs <- boundary_preds[k + 1, idx] - boundary_preds[k, idx]
      avg <- mean(diffs)
    }
    ale_vals[k + 1] <- ale_vals[k] + avg
  }
  bin_id <- findInterval(xj, bins, rightmost.closed = TRUE, all.inside = TRUE)
  bin_id[bin_id == 0] <- 1
  bin_id[bin_id >= length(bins)] <- length(bins) - 1
  ale_at_data <- numeric(length(xj))
  for (i in seq_along(xj)) {
    k <- bin_id[i]
    lower <- bins[k]
    upper <- bins[k + 1]
    if (upper - lower < 1e-8) {
      ale_at_data[i] <- ale_vals[k]
    } else {
      w <- (xj[i] - lower)/(upper - lower)
      ale_at_data[i] <- ale_vals[k] + w * (ale_vals[k + 1] - ale_vals[k])
    }
  }
  ale_vals <- ale_vals - mean(ale_at_data)
  mids <- (bins[-1] + bins[-length(bins)]) / 2
  ale_curve <- (ale_vals[-1] + ale_vals[-length(ale_vals)]) / 2
  data.frame(x = mids, ale = ale_curve)
}

main <- function() {
  root <- resolve_root()
  options(sebart.root = root)
  source(file.path(root, "sebart.R"))
  source(file.path(root, "R", "experiment_utils.R"))

  defaults <- list(
    scenarios = "sigmoid",
    n_train = "200",
    n_test = "1000",
    p = "10",
    bins = "20",
    seed = "202401",
    outdir = file.path(root, "results", "simulations"),
    prefix = "sim_ale",
    standardize = "FALSE"
  )

  args <- modifyList(defaults, parse_args(commandArgs(trailingOnly = TRUE)))
  scenario_vec <- trimws(strsplit(args$scenarios, ",")[[1]])
  n_train <- as.integer(args$n_train)
  n_test <- as.integer(args$n_test)
  p <- as.integer(args$p)
  bins <- max(4L, as.integer(args$bins))
  seed <- as.integer(args$seed)
  outdir <- normalizePath(args$outdir, mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prefix <- args$prefix %||% "sim_ale"
  standardize <- as_bool(args$standardize)

  ale_records <- list()

  for (idx in seq_along(scenario_vec)) {
    scenario <- scenario_vec[idx]
    data <- simulate_sebart_data(
      n_train = n_train,
      n_test = n_test,
      p = p,
      scenario = scenario,
      seed = seed + idx - 1
    )
    X <- as.matrix(data$X_train)
    y <- as.numeric(data$y_train)
    transform_info <- get_transformation(scenario)
    feature_idx <- seq_len(ncol(X))

    design <- build_ale_design(X, feature_idx, bins)
    fit <- sebart(
      y = y,
      X = X,
      X_test = design$X_test,
      ntree = 200,
      nsave = 1000,
      nburn = 1000,
      verbose = FALSE
    )
    sebart_pred <- colMeans(fit$post_ypred)

    bart_fit <- dbarts::bart(
      x.train = X,
      y.train = y,
      x.test = design$X_test,
      ntree = 200,
      ndpost = 1000,
      nskip = 1000,
      verbose = FALSE
    )
    bart_pred <- colMeans(bart_fit$yhat.test)
    oracle_pred <- oracle_observed_mean(
      design$X_test,
      transform_info,
      seed = seed + idx - 1
    )

    block_preds <- lapply(names(design$block_rows), function(name) {
      idxs <- design$block_rows[[name]]
      list(
        sebart = sebart_pred[idxs],
        bart = bart_pred[idxs],
        oracle = oracle_pred[idxs]
      )
    })
    names(block_preds) <- names(design$block_rows)

    for (j in feature_idx) {
      info <- design$feature_info[[as.character(j)]]
      if (is.null(info)) next
      bins_j <- info$bins
      block_names <- info$blocks
      for (model_name in c("sebart", "bart", "oracle")) {
        boundary_preds <- do.call(rbind, lapply(block_names, function(bn) block_preds[[bn]][[model_name]]))
        ale_df <- compute_ale(X[, j], bins_j, boundary_preds)
        if (is.null(ale_df)) next
        ale_df$Scenario <- tools::toTitleCase(scenario)
        ale_df$Feature <- paste0("X", j)
        ale_df$Model <- model_name
        ale_records[[length(ale_records) + 1]] <- ale_df
      }
    }
  }

  if (!length(ale_records)) {
    message("No ALE curves generated.")
    quit(save = "no", status = 0)
  }

  ale_data <- do.call(rbind, ale_records)
  ale_data$Feature <- factor(ale_data$Feature, levels = unique(ale_data$Feature))

  if (standardize) {
    # Standardize ALE values per Scenario x Feature x Model so curves are
    # directly comparable by shape (location/scale removed).
    group_key <- interaction(ale_data$Scenario, ale_data$Feature, ale_data$Model, drop = TRUE)
    ale_data$ale <- ave(
      ale_data$ale,
      group_key,
      FUN = function(v) {
        s <- stats::sd(v, na.rm = TRUE)
        if (!is.finite(s) || s < 1e-12) return(v - mean(v, na.rm = TRUE))
        (v - mean(v, na.rm = TRUE)) / s
      }
    )
  }

  ale_csv <- file.path(outdir, paste0(prefix, "_curves.csv"))
  utils::write.csv(ale_data, ale_csv, row.names = FALSE)

  scenario_label <- paste(unique(ale_data$Scenario), collapse = ", ")
  ale_plot <- ggplot(ale_data, aes(x = x, y = ale, colour = Model)) +
    geom_line(linewidth = 0.9) +
    facet_grid(Scenario ~ Feature, scales = "free") +
    labs(title = sprintf("Accumulated local effects (%s)", scenario_label),
         x = "Feature value",
         y = if (standardize) "Standardized ALE (z-score)" else "ALE on observed scale",
         colour = "Model") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text.x = element_text(size = 8))

  ale_png <- file.path(outdir, paste0(prefix, "_curves.png"))
  ggsave(ale_png, ale_plot, width = 9, height = 6, dpi = 300, bg = "white")

  message("Saved ALE data to ", ale_csv)
  message("Saved ALE plot to ", ale_png)
}

if (sys.nframe() == 0) {
  main()
}
