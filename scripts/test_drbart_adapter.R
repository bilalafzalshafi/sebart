#!/usr/bin/env Rscript
# Lightweight smoke test to verify the DR-BART adapter can fit, predict, and
# return standardized outputs on a toy dataset with multiple covariates.

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

simulate_readme_data <- function(n_train = 200, n_test = 32, seed = 1234) {
  set.seed(seed)
  gamma_shape <- function(x) 0.5 + x ^ 2
  m_fun <- function(x) 1 + 2 * (x - 0.8)
  p_fun <- function(x) exp(-10 * (x - 0.8) ^ 2)
  mu0 <- function(x) 5 * exp(15 * (x - 0.5)) / (1 + exp(15 * (x - 0.5))) - 4 * x
  r_fun <- function(x) {
    n <- length(x)
    z <- rbinom(n, 1, p_fun(x))
    z * rnorm(n, m_fun(x), 0.3) + (1 - z) * log(rgamma(n, gamma_shape(x), 1)) + mu0(x)
  }

  x_train <- runif(n_train)
  x_test <- runif(n_test)
  aux_train <- rnorm(n_train)
  aux_test <- rnorm(n_test)

  list(
    X_train = cbind(x_train, aux_train),
    X_test = cbind(x_test, aux_test),
    y_train = r_fun(x_train)
  )
}

main <- function() {
  root <- resolve_root()
  options(sebart.root = root)
  source(file.path(root, "R", "drbart_adapter.R"))

  if (!requireNamespace("drbart", quietly = TRUE)) {
    stop("Package 'drbart' is not installed. Run scripts/install_drbart.R first.", call. = FALSE)
  }

  data <- simulate_readme_data(n_train = 150, n_test = 25, seed = 2025)
  fit <- drbart_adapter(
    y_train = data$y_train,
    X_train = data$X_train,
    X_test = data$X_test,
    options = list(
      mcmc = list(nburn = 500, nsim = 200, nthin = 1, variance = "ux"),
      ygrid_length = 200,
      chunk_size = 8,
      interval_levels = c(0.8, 0.9),
      predict = list(n_cores = NULL)
    )
  )

  stopifnot(length(fit$fitted.values) == nrow(data$X_test))
  stopifnot(is.null(fit$post_ypred))
  stopifnot(length(fit$quantile_probs) == dim(fit$quantile_predictions)[2])
  stopifnot(dim(fit$quantile_predictions)[1] == nrow(data$X_test))
  stopifnot(length(fit$quantile_probs) > 0)

  cat("Multivariate DR-BART adapter smoke test completed successfully.\n")
}

if (sys.nframe() == 0) {
  main()
}
