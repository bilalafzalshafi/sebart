# Helper utilities to interface with vittorioorlandi/drbart without touching
# internal sampler objects. The adapter keeps the surface area small so the rest
# of the study code can treat DR-BART as another model that returns point
# predictions, posterior predictive draws, and (optionally) densities on a grid.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

default_drbart_adapter_options <- function() {
  list(
    mcmc = list(
      nburn = 1000L,
      nsim = 1000L,
      nthin = 1L,
      m_mean = 200L,
      m_var = 100L,
      alpha = 0.95,
      beta = 2,
      lambda = 1,
      nu = 2,
      kfac = 2,
      phi0 = 1,
      variance = "ux",
      printevery = NULL
    ),
    predict = list(n_cores = NULL),
    ygrid = NULL,
    ygrid_limits = NULL,
    ygrid_length = 256L,
    ygrid_expand = 0.1,
    chunk_size = 64L,
    interval_levels = c(0.8, 0.9, 0.95),
    sample_seed = NULL
  )
}

coerce_predictor_matrix <- function(x, ncol_target = NULL) {
  if (is.null(x)) {
    stop("Predictor matrix cannot be NULL", call. = FALSE)
  }
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("Predictor matrix must be numeric", call. = FALSE)
  }
  if (!is.null(ncol_target) && ncol(x) != ncol_target) {
    stop("Predictor matrices must have the same number of columns", call. = FALSE)
  }
  if (ncol(x) == 1) {
    colnames(x) <- colnames(x) %||% "x1"
  }
  x
}

resolve_y_grid <- function(y, grid = NULL, limits = NULL,
                           length_out = 256L, expand = 0.1) {
  if (!is.null(grid)) {
    grid <- sort(unique(as.numeric(grid)))
    return(grid)
  }
  if (!is.null(limits)) {
    limits <- as.numeric(limits)
    if (length(limits) != 2L) stop("ygrid_limits must have length 2", call. = FALSE)
    ymin <- limits[1]
    ymax <- limits[2]
  } else {
    ymin <- min(y, na.rm = TRUE)
    ymax <- max(y, na.rm = TRUE)
    if (!is.finite(ymin) || !is.finite(ymax)) {
      stop("Cannot determine y-grid range: non-finite values in response", call. = FALSE)
    }
    if (ymin == ymax) {
      buffer <- ifelse(ymin == 0, 0.1, abs(ymin) * expand)
      ymin <- ymin - buffer
      ymax <- ymax + buffer
    } else {
      span <- ymax - ymin
      buffer <- span * expand
      ymin <- ymin - buffer
      ymax <- ymax + buffer
    }
  }
  seq(from = ymin, to = ymax, length.out = length_out)
}

drbart_adapter <- function(y_train, X_train, X_test, options = list()) {
  if (!requireNamespace("drbart", quietly = TRUE)) {
    stop("Package 'drbart' is required. Run scripts/install_drbart.R first.", call. = FALSE)
  }
  opts <- modifyList(default_drbart_adapter_options(), options)
  X_train <- coerce_predictor_matrix(X_train)
  X_test <- coerce_predictor_matrix(X_test, ncol_target = ncol(X_train))
  y_train <- as.numeric(y_train)
  rng_state <- NULL
  if (!is.null(opts$sample_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rng_state <- get(".Random.seed", envir = .GlobalEnv)
    }
    set.seed(opts$sample_seed)
    on.exit({
      if (is.null(rng_state)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", rng_state, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  ygrid <- resolve_y_grid(
    y = y_train,
    grid = opts$ygrid,
    limits = opts$ygrid_limits,
    length_out = opts$ygrid_length,
    expand = opts$ygrid_expand
  )

  fit_args <- modifyList(default_drbart_adapter_options()$mcmc, opts$mcmc)
  variance <- match.arg(fit_args$variance, c("ux", "x", "const"))
  mean_file <- fit_args$mean_file %||% tempfile("drbart_mean_", fileext = ".txt")
  prec_file <- fit_args$prec_file %||% tempfile("drbart_prec_", fileext = ".txt")
  fit_args$variance <- variance
  fit_args$mean_file <- mean_file
  if (variance != "const") {
    fit_args$prec_file <- prec_file
  } else {
    fit_args$prec_file <- NULL
  }
  fit_args <- fit_args[!vapply(fit_args, is.null, logical(1))]

  on.exit({
    files_to_remove <- c(mean_file, prec_file)
    files_to_remove <- files_to_remove[file.exists(files_to_remove)]
    unlink(files_to_remove)
  }, add = TRUE)

  fit <- do.call(drbart::drbart, c(list(y = y_train, x = X_train), fit_args))
  nsim <- length(fit$fit$ucuts)
  if (nsim <= 0) {
    stop("DR-BART fit did not return any posterior samples", call. = FALSE)
  }
  n_test <- nrow(X_test)

  mean_args <- list(object = fit, xpred = X_test, ygrid = ygrid, type = "mean")
  if (!is.null(opts$predict$n_cores)) {
    mean_args$n_cores <- opts$predict$n_cores
  }
  mean_pred <- do.call(stats::predict, mean_args)
  fitted_values <- apply(mean_pred$preds, 1, mean, na.rm = TRUE)

  levels <- sort(unique(opts$interval_levels))
  quant_probs <- sort(unique(c((1 - levels) / 2, 1 - (1 - levels) / 2)))
  quant_args <- list(
    object = fit,
    xpred = X_test,
    ygrid = ygrid,
    type = "quantiles",
    quantiles = quant_probs
  )
  if (!is.null(opts$predict$n_cores)) {
    quant_args$n_cores <- opts$predict$n_cores
  }
  quant_pred <- do.call(stats::predict, quant_args)
  quantile_array <- quant_pred$preds
  interval_summary <- NULL

  list(
    fitted.values = fitted_values,
    post_ypred = NULL,
    y_grid = ygrid,
    quantile_probs = quant_probs,
    quantile_predictions = quantile_array,
    interval_summary = interval_summary,
    model = "drbart"
  )
}
