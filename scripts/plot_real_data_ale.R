#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dbarts)
})

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

parse_dataset_targets <- function(spec) {
  entries <- strsplit(spec, ",")[[1]]
  lapply(entries, function(entry) {
    parts <- strsplit(entry, ":", fixed = TRUE)[[1]]
    dataset <- parts[1]
    target <- if (length(parts) > 1) parts[2] else NA_character_
    list(dataset = dataset, target = target)
  })
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)) a else b

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
  root_dir <- resolve_root()
  options(sebart.root = root_dir)
  source(file.path(root_dir, "R", "experiment_utils.R"))
  source(file.path(root_dir, "sebart.R"))

  defaults <- list(
    datasets = "forestfires:area",
    bins = "20",
    outdir = file.path(root_dir, "results", "real_data"),
    prefix = "real_ale",
    nsave = "500",
    nburn = "500"
  )
  args <- modifyList(defaults, parse_args(commandArgs(trailingOnly = TRUE)))
  dataset_targets <- parse_dataset_targets(args$datasets)
  bins <- max(4L, as.integer(args$bins))
  outdir <- normalizePath(args$outdir, mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  prefix <- args$prefix %||% "real_ale"
  nsave <- as.integer(args$nsave)
  nburn <- as.integer(args$nburn)

  ale_records <- list()

  for (dt in dataset_targets) {
    dataset_name <- dt$dataset
    target_feature <- if (!is.na(dt$target)) dt$target else NULL
    data <- load_real_dataset(dataset_name, target_feature, scale_y = FALSE, scale_X = FALSE)
    X <- as.matrix(data$X)
    y <- as.numeric(data$y)
    feature_idx <- seq_len(ncol(X))

    design <- build_ale_design(X, feature_idx, bins)
    fit <- sebart(
      y = y,
      X = X,
      X_test = design$X_test,
      ntree = 200,
      nsave = nsave,
      nburn = nburn,
      verbose = FALSE
    )
    sebart_pred <- colMeans(fit$post_ypred)

    bart_fit <- dbarts::bart(
      x.train = X,
      y.train = y,
      x.test = design$X_test,
      ntree = 200,
      ndpost = nsave,
      nskip = nburn,
      verbose = FALSE
    )
    bart_pred <- colMeans(bart_fit$yhat.test)

    block_preds <- lapply(names(design$block_rows), function(name) {
      idx <- design$block_rows[[name]]
      list(
        sebart = sebart_pred[idx],
        bart = bart_pred[idx]
      )
    })
    names(block_preds) <- names(design$block_rows)

    for (j in feature_idx) {
      info <- design$feature_info[[as.character(j)]]
      if (is.null(info)) next
      bins_j <- info$bins
      block_names <- info$blocks
      for (model_name in c("sebart", "bart")) {
        boundary_preds <- do.call(rbind, lapply(block_names, function(bn) block_preds[[bn]][[model_name]]))
        ale_df <- compute_ale(X[, j], bins_j, boundary_preds)
        if (is.null(ale_df)) next
        ale_df$Dataset <- dataset_name
        ale_df$Target <- target_feature %||% "default"
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
  ale_data$Dataset <- tools::toTitleCase(ale_data$Dataset)
  ale_data$Feature <- factor(ale_data$Feature, levels = unique(ale_data$Feature))

  ale_csv <- file.path(outdir, paste0(prefix, "_curves.csv"))
  utils::write.csv(ale_data, ale_csv, row.names = FALSE)

  ale_plot <- ggplot(ale_data, aes(x = x, y = ale, colour = Model)) +
    geom_line(linewidth = 0.9) +
    facet_grid(Dataset ~ Feature, scales = "free") +
    labs(title = "Accumulated local effects (Forest Fires)",
         x = "Feature value",
         y = "ALE on observed scale",
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
