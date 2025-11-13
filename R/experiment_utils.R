# Shared utilities for SeBART simulations and real-data studies

.sebart_get_project_root <- function() {
  opt_root <- getOption("sebart.root")
  if (!is.null(opt_root)) {
    return(normalizePath(opt_root))
  }

  file_path <- NULL
  frames <- sys.frames()
  for (idx in rev(seq_along(frames))) {
    candidate <- frames[[idx]]$ofile
    if (!is.null(candidate)) {
      file_path <- candidate
      break
    }
  }

  if (!is.null(file_path)) {
    dir_path <- dirname(normalizePath(file_path))
    if (basename(dir_path) %in% c("R", "scripts")) {
      return(normalizePath(file.path(dir_path, "..")))
    }
    return(dir_path)
  }

  normalizePath(getwd())
}

sebart_data_path <- function(...) {
  file.path(.sebart_get_project_root(), "data", ...)
}

ensure_packages <- function(pkgs) {
  missing <- vapply(pkgs, function(pkg) !requireNamespace(pkg, quietly = TRUE), logical(1))
  if (any(missing)) {
    missing_pkgs <- pkgs[missing]
    stop(
      sprintf(
        "Missing packages: %s. Install them before running the study.",
        paste(missing_pkgs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

compute_ise <- function(estimate, truth) {
  if (length(estimate) != length(truth)) {
    stop("estimate and truth must have the same length", call. = FALSE)
  }
  mean((estimate - truth)^2)
}

compute_transformation_diagnostics <- function(draws, truth, y_grid,
                                               mu = NULL, sigma = NULL,
                                               record_curves = FALSE,
                                               curve_levels = c(0.05, 0.95)) {
  if (is.null(draws) || is.null(truth)) {
    return(NULL)
  }
  if (ncol(draws) != length(truth)) {
    stop("Transformation draws and truth have incompatible dimensions", call. = FALSE)
  }

  draws_raw <- draws
  truth_raw <- matrix(truth, nrow = nrow(draws_raw), ncol = length(truth), byrow = TRUE)

  raw_mean <- colMeans(draws_raw)
  raw_ise_mean <- compute_ise(raw_mean, truth)
  raw_ise_draw <- mean(rowMeans((draws_raw - truth_raw)^2))

  result <- list(
    ise_raw = raw_ise_mean,
    ise_raw_draw = raw_ise_draw,
    mean_raw = raw_mean,
    mean_truth_raw = truth
  )

  if (record_curves) {
    lower_prob <- min(curve_levels)
    upper_prob <- max(curve_levels)
    raw_lower <- apply(draws_raw, 2, quantile, probs = lower_prob, na.rm = TRUE)
    raw_upper <- apply(draws_raw, 2, quantile, probs = upper_prob, na.rm = TRUE)

    result$curves <- rbind(
      data.frame(
        y = y_grid,
        truth = truth,
        mean = raw_mean,
        lower = raw_lower,
        upper = raw_upper,
        stringsAsFactors = FALSE
      )
    )
  }

  result
}

# -------------------- Simulation utilities --------------------

simulate_sebart_data <- function(n_train = 200, n_test = 1000, p = 10,
                                 scenario = "sigmoid", seed = NULL,
                                 predictor = c("gaussian", "uniform")) {
  predictor <- match.arg(predictor)
  if (!is.null(seed)) set.seed(seed)

  draw_predictors <- switch(
    predictor,
    gaussian = function(n) matrix(rnorm(n * p), n, p),
    uniform = function(n) matrix(runif(n * p), n, p)
  )

  X_train <- draw_predictors(n_train)
  X_test <- draw_predictors(n_test)

  f_true <- function(X) {
    10 * sin(pi * X[, 1] * X[, 2]) +
      20 * (X[, 3] - 0.5)^2 +
      10 * X[, 4] + 5 * X[, 5]
  }

  f_true_train <- f_true(X_train)
  f_true_test <- f_true(X_test)

  z_train <- f_true_train + rnorm(n_train, sd = 1)
  z_test <- f_true_test + rnorm(n_test, sd = 1)

  transform_info <- get_transformation(scenario)
  y_train <- transform_info$g_inv(z_train)
  y_test <- transform_info$g_inv(z_test)

  y_unique <- sort(unique(y_train))
  g_true <- if (is.null(transform_info$g)) NULL else transform_info$g(y_unique)

  list(
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test,
    f_true_train = f_true_train,
    f_true_test = f_true_test,
    z_train = z_train,
    z_test = z_test,
    scenario = scenario,
    g_inverse = transform_info$g_inv,
    g_true = g_true,
    y_unique = y_unique,
    description = transform_info$description
  )
}

get_transformation <- function(scenario) {
  transformations <- list(
    step = list(
      g = function(y) y,
      g_inv = function(z) {
        g_steps <- rexp(n = 10)
        approxfun(
          seq(-3, 3, length.out = 10),
          cumsum(g_steps),
          rule = 2
        )((z - mean(z)) / sd(z))
      },
      description = "Monotonic step transformation"
    ),
    sigmoid = list(
      g = function(y) -2 * log(20 / pmax(y, 1e-10) - 1),
      g_inv = function(z) 20 / (1 + exp(-z / 2)),
      description = "Inverse sigmoid transformation, bounded [0,20]"
    ),
    beta = list(
      g = function(y) qnorm(pbeta(pmax(pmin(y, 1 - 1e-10), 1e-10), 0.1, 0.5)),
      g_inv = function(z) {
        z_std <- (z - mean(z)) / sd(z)
        qbeta(pnorm(z_std), shape1 = 0.1, shape2 = 0.5)
      },
      description = "Beta(0.1,0.5) bounded [0,1], many values near zero"
    ),
    arctangent = list(
      g = function(y) 5 * tan(y / 100),
      g_inv = function(z) 100 * atan(z / 5),
      description = "Arctangent transformation, smooth bounded"
    ),
    box_cox = list(
      g = function(y) SeBR:::g_bc(y, lambda = 0.5),
      g_inv = function(z) SeBR:::g_inv_bc(z, lambda = 0.5),
      description = "Inverse Box-Cox square root transformation"
    ),
    almost_linear = list(
      g = function(y) (-1 + sqrt(1 + 0.4 * y)) / 0.2,
      g_inv = function(z) z + 0.1 * z^2,
      description = "Nearly linear with weak quadratic component"
    ),
    polynomial = list(
      g = function(y) y^(1 / 2.5),
      g_inv = function(z) {
        z_positive <- z - min(z) + 5
        z_positive^2.5
      },
      description = "Polynomial x^2.5 transformation"
    )
  )

  if (!scenario %in% names(transformations)) {
    stop(
      "Unknown scenario. Available:",
      paste(names(transformations), collapse = ", "),
      call. = FALSE
    )
  }

  transformations[[scenario]]
}

get_available_scenarios <- function() {
  c("sigmoid", "beta", "arctangent", "box_cox", "almost_linear", "polynomial", "step")
}

fit_models <- function(y_train, X_train, X_test, models = c("sebart", "bart", "bart_bc", "sblm"),
                       config = list(ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE),
                       progress = NULL) {
  ensure_packages("SeBR")
  fits <- list()
  total_models <- length(models)

  emit_progress <- function(model_idx, model_name) {
    if (is.null(progress) || !isTRUE(progress$verbose)) return()
    iter_txt <- NULL
    if (!is.null(progress$iteration) && !is.null(progress$total_iterations)) {
      iter_txt <- sprintf("Iteration %d/%d", progress$iteration, progress$total_iterations)
    }
    model_txt <- sprintf("Model %d/%d: %s", model_idx, total_models, model_name)
    if (!is.null(iter_txt)) {
      message(sprintf("%s, %s", iter_txt, model_txt))
    } else {
      message(model_txt)
    }
  }

  for (model_idx in seq_along(models)) {
    model_name <- models[model_idx]
    emit_progress(model_idx, model_name)
    fit <- tryCatch(
      {
        runtime <- system.time({
          result <- switch(
            model_name,
            sebart = {
              ensure_packages("dbarts")
              sebart(
                y = y_train,
                X = X_train,
                X_test = X_test,
                ntree = config$ntree,
                nsave = config$nsave,
                nburn = config$nburn,
                verbose = config$verbose
              )
            },
            bart = {
              ensure_packages("dbarts")
              fit <- dbarts::bart(
                x.train = X_train,
                y.train = y_train,
                x.test = X_test,
                ntree = config$ntree,
                ndpost = config$nsave,
                nskip = config$nburn,
                verbose = config$verbose
              )
              list(
                fitted.values = fit$yhat.test.mean,
                post_ypred = fit$yhat.test,
                model = "bart"
              )
            },
            bart_bc = {
              bart_bc(
                y = y_train,
                X = X_train,
                X_test = X_test,
                ntree = config$ntree,
                nsave = config$nsave,
                nburn = config$nburn,
                verbose = config$verbose
              )
            },
            sblm = {
              SeBR::sblm(
                y = y_train,
                X = X_train,
                X_test = X_test,
                verbose = FALSE
              )
            },
            lm = {
              df_train <- data.frame(y = y_train, X_train)
              colnames(df_train) <- c("y", paste0("X", seq_len(ncol(X_train))))
              lm_fit <- stats::lm(y ~ ., data = df_train)
              preds <- as.numeric(stats::predict(lm_fit, newdata = data.frame(X_test)))
              resid_sd <- sqrt(mean(stats::residuals(lm_fit)^2, na.rm = TRUE))
              draws <- matrix(
                rnorm(config$nsave * length(preds), mean = rep(preds, each = config$nsave), sd = resid_sd),
                nrow = config$nsave,
                byrow = TRUE
              )
              list(
                fitted.values = preds,
                post_ypred = draws,
                model = "lm",
                runtime = NA_real_
              )
            },
            drbart = {
              if (!requireNamespace("drbart", quietly = TRUE)) {
                warning("Skipping drbart: package not installed.")
                NULL
              } else {
                fit_drbart(y_train, X_train, X_test)
              }
            },
            stop("Unknown model: ", model_name)
          )
        })
        if (!is.null(result)) {
          result$runtime <- as.numeric(runtime["elapsed"])
        }
        result
      },
      error = function(e) {
        warning(sprintf("Model %s failed: %s", model_name, e$message))
        NULL
      }
    )
    fits[[model_name]] <- fit
  }
  fits
}

fit_drbart <- function(y, X, X_test,
                       nburn = 1000, nsim = 1000, nthin = 1,
                       m_mean = 200, m_var = 100,
                       variance = c("ux", "x", "const")) {
  variance <- match.arg(variance)
  X <- as.matrix(X)
  X_test <- as.matrix(X_test)

  mean_file <- tempfile("drbart_mean_", fileext = ".txt")
  prec_file <- tempfile("drbart_prec_", fileext = ".txt")
  on.exit(unlink(c(mean_file, prec_file)), add = TRUE)

  fit <- drbart::drbart(
    y = y,
    x = X,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    m_mean = m_mean,
    m_var = m_var,
    variance = variance,
    mean_file = mean_file,
    prec_file = prec_file
  )

  compute_posterior_means <- function(fit_object, new_x) {
    ts_mean <- drbart:::TreeSamples$new()
    ts_mean$load(mean_file)

    fit_internal <- fit_object$fit
    nsim_local <- length(fit_internal$ucuts)
    if (nsim_local == 0) {
      stop("No posterior samples available from drbart fit")
    }

    logprobs <- lapply(fit_internal$ucuts, function(u) log(diff(c(0, u, 1))))
    mids <- lapply(fit_internal$ucuts, function(u) c(0, u) + diff(c(0, u, 1)) / 2)

    n_test <- nrow(new_x)
    posterior_means <- matrix(NA_real_, n_test, nsim_local)

    for (j in seq_len(n_test)) {
      row_values <- as.list(new_x[j, ])
      des <- lapply(mids, function(m) do.call(data.frame, c(list(m), row_values)))

      for (i in seq_len(nsim_local)) {
        mu_vec <- c(ts_mean$predict_i(t(des[[i]]), i - 1))
        weights <- exp(logprobs[[i]])
        posterior_means[j, i] <- sum(weights * mu_vec)
      }
    }

    posterior_means
  }

  posterior_means <- tryCatch(
    compute_posterior_means(fit, X_test),
    error = function(e) {
      warning("drbart prediction failed: ", conditionMessage(e))
      matrix(NA_real_, nrow(X_test), 1)
    }
  )

  fitted_mean <- rowMeans(posterior_means, na.rm = TRUE)
  fitted_mean[is.nan(fitted_mean)] <- NA_real_

  list(
    fitted.values = fitted_mean,
    post_ypred = NULL,
    model = "drbart"
  )
}

rmse <- function(pred, truth) {
  sqrt(mean((pred - truth)^2, na.rm = TRUE))
}

calc_interval_stats <- function(draws, y_true, alpha = 0.1) {
  if (is.null(draws)) return(list(coverage = NA_real_, width = NA_real_))
  lower <- apply(draws, 2, quantile, alpha / 2, na.rm = TRUE)
  upper <- apply(draws, 2, quantile, 1 - alpha / 2, na.rm = TRUE)
  coverage <- mean(y_true >= lower & y_true <= upper, na.rm = TRUE)
  width <- mean(upper - lower, na.rm = TRUE)
  list(coverage = coverage, width = width)
}

calc_interval_stats_levels <- function(draws, y_true, levels = c(0.8, 0.9, 0.95)) {
  if (is.null(draws)) {
    return(data.frame(
      Level = levels,
      Coverage = NA_real_,
      Width = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  res <- lapply(levels, function(level) {
    alpha <- 1 - level
    lower <- apply(draws, 2, quantile, alpha / 2, na.rm = TRUE)
    upper <- apply(draws, 2, quantile, 1 - alpha / 2, na.rm = TRUE)
    data.frame(
      Level = level,
      Coverage = mean(y_true >= lower & y_true <= upper, na.rm = TRUE),
      Width = mean(upper - lower, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

calc_hpd_stats_levels <- function(draws, y_true, levels = c(0.8, 0.9, 0.95)) {
  if (is.null(draws)) {
    return(data.frame(
      Level = levels,
      Coverage = NA_real_,
      Width = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  ensure_packages("coda")
  res <- lapply(levels, function(level) {
    alpha <- 1 - level
    intervals <- t(vapply(seq_len(ncol(draws)), function(j) {
      sample_vec <- draws[, j]
      if (all(is.na(sample_vec))) {
        c(NA_real_, NA_real_)
      } else {
        hpd <- coda::HPDinterval(coda::mcmc(sample_vec), prob = level)
        c(hpd[1], hpd[2])
      }
    }, numeric(2)))
    lower <- intervals[, 1]
    upper <- intervals[, 2]
    data.frame(
      Level = level,
      Coverage = mean(y_true >= lower & y_true <= upper, na.rm = TRUE),
      Width = mean(upper - lower, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

crps_empirical <- function(samples, obs) {
  if (length(samples) == 0L || all(is.na(samples))) return(NA_real_)
  samples <- samples[is.finite(samples)]
  if (length(samples) == 0L) return(NA_real_)
  term1 <- mean(abs(samples - obs))
  sx <- sort(samples)
  n <- length(sx)
  coeff <- (2 * seq_len(n) - n - 1)
  term2 <- sum(coeff * sx) / (n^2)
  term1 - term2
}

compute_crps_vector <- function(draws, y_true) {
  if (is.null(draws)) return(rep(NA_real_, length(y_true)))
  apply(draws, 2, crps_empirical, obs = y_true)
}

sample_posterior_predictions <- function(draws, sample_size = 500) {
  if (is.null(draws)) return(NULL)
  values <- as.numeric(draws)
  values <- values[is.finite(values)]
  if (!length(values)) return(NULL)
  sample(values, size = min(sample_size, length(values)), replace = length(values) < sample_size)
}

evaluate_model <- function(fit_object, y_test, model_name, alpha = 0.1) {
  if (is.null(fit_object)) {
    return(data.frame(
      Model = model_name,
      RMSE = NA_real_,
      Coverage_Q_90 = NA_real_,
      Width_Q_90 = NA_real_,
      Coverage_Q_80 = NA_real_,
      Width_Q_80 = NA_real_,
      Coverage_Q_95 = NA_real_,
      Width_Q_95 = NA_real_,
      Coverage_HPD_80 = NA_real_,
      Width_HPD_80 = NA_real_,
      Coverage_HPD_90 = NA_real_,
      Width_HPD_90 = NA_real_,
      Coverage_HPD_95 = NA_real_,
      Width_HPD_95 = NA_real_,
      Mean_CRPS = NA_real_,
      Median_CRPS = NA_real_,
      Runtime = NA_real_
    ))
  }

  preds <- fit_object$fitted.values
  rmse_val <- rmse(preds, y_test)
  levels <- c(0.8, 0.9, 0.95)
  level_labels <- sprintf("%02d", levels * 100)

  interval_q <- calc_interval_stats_levels(fit_object$post_ypred, y_test, levels)
  interval_hpd <- calc_hpd_stats_levels(fit_object$post_ypred, y_test, levels)
  crps_vals <- compute_crps_vector(fit_object$post_ypred, y_test)

  df <- data.frame(
    Model = model_name,
    RMSE = rmse_val,
    Mean_CRPS = mean(crps_vals, na.rm = TRUE),
    Median_CRPS = median(crps_vals, na.rm = TRUE),
    Runtime = fit_object$runtime,
    stringsAsFactors = FALSE
  )

  for (idx in seq_along(levels)) {
    label <- level_labels[idx]
    df[[paste0("Coverage_Q_", label)]] <- interval_q$Coverage[idx]
    df[[paste0("Width_Q_", label)]] <- interval_q$Width[idx]
    df[[paste0("Coverage_HPD_", label)]] <- interval_hpd$Coverage[idx]
    df[[paste0("Width_HPD_", label)]] <- interval_hpd$Width[idx]
  }

  df
}

run_simulation_iteration <- function(scenario, iteration, total_iterations, config, models,
                                     n_train, n_test, p, seed_offset = 0,
                                     alpha = 0.1,
                                     record_transformation = FALSE,
                                     verbose_progress = FALSE) {
  seed <- config$seed + iteration - 1 + seed_offset
  data <- simulate_sebart_data(
    n_train = n_train,
    n_test = n_test,
    p = p,
    scenario = scenario,
    seed = seed,
    predictor = config$predictor
  )

  fits <- fit_models(
    y_train = data$y_train,
    X_train = data$X_train,
    X_test = data$X_test,
    models = models,
    config = config$model,
    progress = list(
      verbose = verbose_progress,
      iteration = iteration,
      total_iterations = total_iterations
    )
  )

  metrics <- do.call(rbind, lapply(names(fits), function(model_name) {
    evaluate_model(fits[[model_name]], data$y_test, model_name, alpha)
  }))
  metrics$Scenario <- scenario
  metrics$Iteration <- iteration
  metrics$N_train <- n_train

  action_range <- data.frame(
    Scenario = scenario,
    Iteration = iteration,
    N_train = n_train,
    Latent_Min = min(c(data$z_train, data$z_test)),
    Latent_Max = max(c(data$z_train, data$z_test)),
    stringsAsFactors = FALSE
  )

  transformation_metrics <- NULL
  transformation_curves <- NULL
  if (!is.null(data$g_true) && "sebart" %in% names(fits) && !is.null(fits$sebart$post_g_raw)) {
    diag <- compute_transformation_diagnostics(
      draws = fits$sebart$post_g_raw,
      truth = data$g_true,
      y_grid = data$y_unique,
      mu = fits$sebart$post_g_mu,
      sigma = fits$sebart$post_g_sigma,
      record_curves = record_transformation
    )
    if (!is.null(diag)) {
      transformation_metrics <- data.frame(
        Scenario = scenario,
        Iteration = iteration,
        N_train = n_train,
        ISE_raw_mean = diag$ise_raw,
        ISE_raw_expected = diag$ise_raw_draw,
        Mu_mean = if (!is.null(fits$sebart$post_g_mu)) mean(fits$sebart$post_g_mu) else NA_real_,
        Sigma_mean = if (!is.null(fits$sebart$post_g_sigma)) mean(fits$sebart$post_g_sigma) else NA_real_,
        stringsAsFactors = FALSE
      )
      if (!is.null(diag$curves)) {
        transformation_curves <- cbind(
          Scenario = scenario,
          Iteration = iteration,
          N_train = n_train,
          diag$curves
        )
      }
    }
  }

  list(
    metrics = metrics,
    action_range = action_range,
    transformation_metrics = transformation_metrics,
    transformation_curves = transformation_curves
  )
}

summarize_simulation_metrics <- function(metrics_df) {
  dplyr::summarise(
    dplyr::group_by(metrics_df, Scenario, Model, N_train),
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    Coverage_Q_80 = mean(Coverage_Q_80, na.rm = TRUE),
    Width_Q_80 = mean(Width_Q_80, na.rm = TRUE),
    Coverage_Q_90 = mean(Coverage_Q_90, na.rm = TRUE),
    Width_Q_90 = mean(Width_Q_90, na.rm = TRUE),
    Coverage_Q_95 = mean(Coverage_Q_95, na.rm = TRUE),
    Width_Q_95 = mean(Width_Q_95, na.rm = TRUE),
    Coverage_HPD_80 = mean(Coverage_HPD_80, na.rm = TRUE),
    Width_HPD_80 = mean(Width_HPD_80, na.rm = TRUE),
    Coverage_HPD_90 = mean(Coverage_HPD_90, na.rm = TRUE),
    Width_HPD_90 = mean(Width_HPD_90, na.rm = TRUE),
    Coverage_HPD_95 = mean(Coverage_HPD_95, na.rm = TRUE),
    Width_HPD_95 = mean(Width_HPD_95, na.rm = TRUE),
    Mean_CRPS = mean(Mean_CRPS, na.rm = TRUE),
    Median_CRPS = median(Median_CRPS, na.rm = TRUE),
    Mean_Runtime = mean(Runtime, na.rm = TRUE),
    .groups = "drop"
  )
}

plot_interval_widths <- function(metrics_df, scenario_filter = NULL, title = NULL) {
  ensure_packages(c("ggplot2", "dplyr"))
  df <- metrics_df
  if (!is.null(scenario_filter)) {
    df <- df[df$Scenario %in% scenario_filter, , drop = FALSE]
  }
  df <- df[!is.na(df$Width_Q_90), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  coverage_summary <- dplyr::summarise(
    dplyr::group_by(df, Model, Scenario, N_train),
    Mean_Coverage = mean(Coverage_Q_90, na.rm = TRUE),
    .groups = "drop"
  )

  max_width <- max(df$Width_Q_90) * 1.15

  df$Scenario_Label <- paste0(df$Scenario, " (n=", df$N_train, ")")
  coverage_summary$Scenario_Label <- paste0(coverage_summary$Scenario, " (n=", coverage_summary$N_train, ")")

  ggplot2::ggplot(df, ggplot2::aes(x = Width_Q_90, y = Model)) +
    ggplot2::geom_boxplot(fill = "white", colour = "black", outlier.size = 1,
                          outlier.shape = 1, size = 0.5) +
    ggplot2::geom_text(
      data = coverage_summary,
      ggplot2::aes(x = max_width, y = Model,
                   label = paste0(round(Mean_Coverage * 100), "%")),
      colour = "blue", hjust = 0, size = 3.3
    ) +
    ggplot2::facet_wrap(~Scenario_Label, ncol = 1) +
    ggplot2::labs(x = "Interval width", y = NULL,
                  title = title %||% "90% prediction interval widths") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.20)))
}

save_interval_widths_by_scenario <- function(metrics_df, outdir, prefix,
                                             scenarios = NULL,
                                             n_train = NULL,
                                             levels = c("80", "90", "95"),
                                             interval_types = c("Q", "HPD"),
                                             allowed_models = c("sebart", "bart", "bart_bc", "sblm")) {
  ensure_packages(c("ggplot2", "dplyr", "tidyr"))
  if (is.null(metrics_df) || !nrow(metrics_df)) return(character())

  df <- metrics_df
  if (!is.null(n_train)) {
    df <- df[df$N_train %in% n_train, , drop = FALSE]
  }
  if (!is.null(scenarios)) {
    df <- df[df$Scenario %in% scenarios, , drop = FALSE]
  }
  df <- df[df$Model %in% allowed_models & !is.na(df$Model) & nzchar(df$Model), , drop = FALSE]
  if (!nrow(df)) return(character())

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  type_levels <- intersect(interval_types, c("Q", "HPD"))
  level_levels <- intersect(levels, c("80", "90", "95"))
  if (!length(type_levels) || !length(level_levels)) return(character())

  df <- df[!is.na(df$Model) & nzchar(df$Model), , drop = FALSE]
  scenario_order <- unique(df$Scenario)
  model_levels <- intersect(allowed_models, unique(df$Model))
  type_labels <- c(Q = "Equal-tail", HPD = "HPD")
  fill_values <- c(
    sebart = "#1d4ed8",
    bart = "#9333ea",
    bart_bc = "#059669",
    sblm = "#f97316"
  )
  fill_values <- fill_values[names(fill_values) %in% model_levels]

  saved <- character()

  for (scenario in scenario_order) {
    scenario_df <- df[df$Scenario == scenario, , drop = FALSE]
    if (!nrow(scenario_df)) next

    width_cols <- grep("^Width_(Q|HPD)_", names(scenario_df), value = TRUE)
    if (!length(width_cols)) next

    long_width <- tidyr::pivot_longer(
      scenario_df,
      cols = dplyr::all_of(width_cols),
      names_to = c("Type", "Level"),
      names_pattern = "Width_(Q|HPD)_(\\d+)",
      values_to = "Width",
      values_drop_na = TRUE
    )

    long_width <- long_width[long_width$Type %in% type_levels & long_width$Level %in% level_levels, , drop = FALSE]
    if (!nrow(long_width)) next

    long_width$Type_Label <- type_labels[long_width$Type]
    long_width$Level_Label <- paste0(long_width$Level, "%")
    long_width$Model <- factor(long_width$Model, levels = model_levels)

    scenario_title <- tools::toTitleCase(scenario)

    coverage_cols <- grep("^Coverage_(Q|HPD)_", names(scenario_df), value = TRUE)
    coverage_long <- tidyr::pivot_longer(
      scenario_df,
      cols = dplyr::all_of(coverage_cols),
      names_to = c("Type", "Level"),
      names_pattern = "Coverage_(Q|HPD)_(\\d+)",
      values_to = "Coverage",
      values_drop_na = TRUE
    )
    coverage_long <- coverage_long[coverage_long$Type %in% type_levels & coverage_long$Level %in% level_levels, , drop = FALSE]
    coverage_long$Type_Label <- type_labels[coverage_long$Type]
    coverage_long$Level_Label <- paste0(coverage_long$Level, "%")
    coverage_long$Model <- factor(coverage_long$Model, levels = model_levels)

    coverage_ann <- dplyr::summarise(
      dplyr::group_by(coverage_long, Model, Type_Label, Level_Label),
      Coverage = mean(Coverage, na.rm = TRUE),
      .groups = "drop"
    )
    width_max <- dplyr::summarise(
      dplyr::group_by(long_width, Model, Type_Label, Level_Label),
      MaxWidth = max(Width, na.rm = TRUE),
      .groups = "drop"
    )
    coverage_ann <- dplyr::left_join(coverage_ann, width_max,
                                     by = c("Model", "Type_Label", "Level_Label"))
    coverage_ann$Label <- sprintf("%s%%", round(coverage_ann$Coverage * 100))
    coverage_ann$MaxWidth <- ifelse(is.finite(coverage_ann$MaxWidth), coverage_ann$MaxWidth, NA_real_)

    plot_obj <- ggplot2::ggplot(long_width, ggplot2::aes(x = Width, y = Model, fill = Model)) +
      ggplot2::geom_boxplot(alpha = 0.35, outlier.size = 0.8, linewidth = 0.3) +
      ggplot2::facet_grid(Type_Label ~ Level_Label, scales = "free_x") +
      ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
      ggplot2::geom_text(
        data = coverage_ann,
        ggplot2::aes(x = MaxWidth, y = Model, label = Label),
        colour = "#1d4ed8",
        hjust = -0.15,
        size = 2.8
      ) +
      ggplot2::labs(
        title = sprintf("Interval widths: %s (n=%s)", scenario_title, paste(unique(scenario_df$N_train), collapse = ",")),
        x = "Interval width",
        y = NULL,
        fill = "Model"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "white", colour = "grey70"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.20)))

    outfile <- file.path(outdir, sprintf("%s_%s_interval_widths.png", prefix, scenario))
    ggplot2::ggsave(outfile, plot_obj, width = 8, height = 4.5, dpi = 300, bg = "white")
    saved <- c(saved, outfile)
  }

  saved
}

save_interval_widths_by_dataset <- function(metrics_df, outdir, prefix,
                                            datasets = NULL,
                                            levels = c("80", "90", "95"),
                                            interval_types = c("Q", "HPD"),
                                            allowed_models = c("sebart", "bart", "bart_bc", "sblm")) {
  ensure_packages(c("ggplot2", "dplyr", "tidyr"))
  if (is.null(metrics_df) || !nrow(metrics_df)) return(character())

  df <- metrics_df
  if (!is.null(datasets)) {
    df <- df[df$Dataset %in% datasets, , drop = FALSE]
  }
  df <- df[df$Model %in% allowed_models & !is.na(df$Model) & nzchar(df$Model), , drop = FALSE]
  if (!nrow(df)) return(character())

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  type_levels <- intersect(interval_types, c("Q", "HPD"))
  level_levels <- intersect(levels, c("80", "90", "95"))
  if (!length(type_levels) || !length(level_levels)) return(character())

  dataset_order <- unique(df$Dataset)
  model_levels <- intersect(allowed_models, unique(df$Model))
  type_labels <- c(Q = "Equal-tail", HPD = "HPD")
  fill_values <- c(
    sebart = "#1d4ed8",
    bart = "#9333ea",
    bart_bc = "#059669",
    sblm = "#f97316"
  )
  fill_values <- fill_values[names(fill_values) %in% model_levels]

  saved <- character()

  for (dataset in dataset_order) {
    data_df <- df[df$Dataset == dataset, , drop = FALSE]
    if (!nrow(data_df)) next

    width_cols = grep("^Width_(Q|HPD)_", names(data_df), value = TRUE)
    if (!length(width_cols)) next

    long_width <- tidyr::pivot_longer(
      data_df,
      cols = dplyr::all_of(width_cols),
      names_to = c("Type", "Level"),
      names_pattern = "Width_(Q|HPD)_(\\d+)",
      values_to = "Width",
      values_drop_na = TRUE
    )

    long_width <- long_width[long_width$Type %in% type_levels & long_width$Level %in% level_levels, , drop = FALSE]
    if (!nrow(long_width)) next

    long_width$Type_Label <- type_labels[long_width$Type]
    long_width$Level_Label <- paste0(long_width$Level, "%")
    long_width$Model <- factor(long_width$Model, levels = model_levels)

    coverage_cols <- grep("^Coverage_(Q|HPD)_", names(data_df), value = TRUE)
    coverage_long <- tidyr::pivot_longer(
      data_df,
      cols = dplyr::all_of(coverage_cols),
      names_to = c("Type", "Level"),
      names_pattern = "Coverage_(Q|HPD)_(\\d+)",
      values_to = "Coverage",
      values_drop_na = TRUE
    )

    coverage_long <- coverage_long[coverage_long$Type %in% type_levels & coverage_long$Level %in% level_levels, , drop = FALSE]
    coverage_long$Type_Label <- type_labels[coverage_long$Type]
    coverage_long$Level_Label <- paste0(coverage_long$Level, "%")
    coverage_long$Model <- factor(coverage_long$Model, levels = model_levels)

    coverage_ann <- dplyr::summarise(
      dplyr::group_by(coverage_long, Model, Type_Label, Level_Label),
      Coverage = mean(Coverage, na.rm = TRUE),
      .groups = "drop"
    )
    width_max <- dplyr::summarise(
      dplyr::group_by(long_width, Model, Type_Label, Level_Label),
      MaxWidth = max(Width, na.rm = TRUE),
      .groups = "drop"
    )
    coverage_ann <- dplyr::left_join(coverage_ann, width_max,
                                     by = c("Model", "Type_Label", "Level_Label"))
    coverage_ann$Label <- sprintf("%s%%", round(coverage_ann$Coverage * 100))
    coverage_ann$MaxWidth <- ifelse(is.finite(coverage_ann$MaxWidth), coverage_ann$MaxWidth, NA_real_)

    plot_obj <- ggplot2::ggplot(long_width, ggplot2::aes(x = Width, y = Model, fill = Model)) +
      ggplot2::geom_boxplot(alpha = 0.35, outlier.size = 0.8, linewidth = 0.3) +
      ggplot2::facet_grid(Type_Label ~ Level_Label, scales = "free_x") +
      ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
      ggplot2::geom_text(
        data = coverage_ann,
        ggplot2::aes(x = MaxWidth, y = Model, label = Label),
        colour = "#1d4ed8",
        hjust = -0.15,
        size = 2.8
      ) +
      ggplot2::labs(
        title = sprintf("Interval widths: %s", tools::toTitleCase(dataset)),
        x = "Interval width",
        y = NULL,
        fill = "Model"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "white", colour = "grey70"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.20)))

    outfile <- file.path(outdir, sprintf("%s_%s_interval_widths.png", prefix, dataset))
    ggplot2::ggsave(outfile, plot_obj, width = 8, height = 4.5, dpi = 300, bg = "white")
    saved <- c(saved, outfile)
  }

  saved
}

prepare_transformation_labels <- function(df) {
  df$Scenario_Label <- sprintf("%s (n=%s)", tools::toTitleCase(df$Scenario), df$N_train)
  df
}

build_truth_panel <- function(df) {
  truth_df <- df[is.finite(df$truth), c("Scenario_Label", "y", "truth")]
  dplyr::distinct(truth_df)
}

make_transformation_plot <- function(df, iteration, prefix) {
  keep <- is.finite(df$y) & is.finite(df$mean) & is.finite(df$lower) & is.finite(df$upper)
  df <- df[keep, , drop = FALSE]
  df <- df[order(df$Scenario, df$y), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  df <- prepare_transformation_labels(df)
  if (!nrow(df)) return(NULL)

  truth_df <- build_truth_panel(df)
  ggplot2::ggplot(df, ggplot2::aes(x = y)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#93c5fd",
      alpha = 0.3,
      colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = mean),
      colour = "#2563eb",
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      data = truth_df,
      ggplot2::aes(y = truth),
      colour = "black",
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~Scenario_Label, scales = "free") +
    ggplot2::labs(
      title = "Transformation recovery",
      subtitle = sprintf("%s — iteration %s", prefix, iteration),
      x = "Observed scale y",
      y = "Estimated g(y)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey70"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
}

save_transformation_plots <- function(curve_df, outdir, prefix,
                                      width = 10, height = 6, dpi = 300) {
  ensure_packages(c("ggplot2", "dplyr"))
  if (is.null(curve_df) || !nrow(curve_df)) return(character())
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  saved_paths <- character()
  iterations <- sort(unique(curve_df$Iteration))
  for (iter in iterations) {
    iter_df <- curve_df[curve_df$Iteration == iter, , drop = FALSE]
    plot_obj <- make_transformation_plot(iter_df, iter, prefix)
    if (is.null(plot_obj)) next
    outfile <- file.path(outdir, sprintf("%s_transformation_iter%s.png", prefix, iter))
    ggplot2::ggsave(outfile, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
    saved_paths <- c(saved_paths, outfile)
  }
  saved_paths
}

run_simulation_study <- function(scenarios,
                                 iterations = 100,
                                 n_train = 200,
                                 n_test = 1000,
                                 p = 10,
                                 models = c("sebart", "bart", "bart_bc", "sblm"),
                                 seed = 202401,
                                 predictor = "gaussian",
                                 alpha = 0.1,
                                 verbose = TRUE,
                                 plot_iterations = 1) {
  ensure_packages(c("dplyr"))
  model_config <- list(ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE)
  config <- list(
    seed = seed,
    predictor = predictor,
    model = list(ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE)
  )

  n_train_vec <- sort(unique(n_train))
  plot_iterations <- min(plot_iterations, iterations)

  metrics_list <- list()
  action_list <- list()
  transformation_metric_list <- list()
  transformation_curve_list <- list()

  for (n_idx in seq_along(n_train_vec)) {
    n_current <- n_train_vec[n_idx]
    for (s_idx in seq_along(scenarios)) {
      scenario <- scenarios[s_idx]
      if (verbose) {
        message(sprintf("Running scenario %s (n_train = %d)", scenario, n_current))
      }

      iter_results <- lapply(seq_len(iterations), function(it) {
        run_simulation_iteration(
          scenario = scenario,
          iteration = it,
          total_iterations = iterations,
          config = config,
          models = models,
          n_train = n_current,
          n_test = n_test,
          p = p,
          alpha = alpha,
          seed_offset = 100000 * (n_idx - 1) + 1000 * (s_idx - 1),
          record_transformation = it <= plot_iterations,
          verbose_progress = verbose
        )
      })

      metrics_list[[paste(scenario, n_current, sep = "_")]] <- dplyr::bind_rows(lapply(iter_results, `[[`, "metrics"))
      action_df <- dplyr::bind_rows(lapply(iter_results, `[[`, "action_range"))
      if (nrow(action_df) > 0) {
        action_list[[paste(scenario, n_current, sep = "_")]] <- action_df
      }

      trans_metric_entries <- lapply(iter_results, `[[`, "transformation_metrics")
      trans_metric_entries <- trans_metric_entries[!vapply(trans_metric_entries, is.null, logical(1))]
      if (length(trans_metric_entries) > 0) {
        transformation_metric_list[[paste(scenario, n_current, sep = "_")]] <- dplyr::bind_rows(trans_metric_entries)
      }

      if (plot_iterations > 0) {
        trans_curves <- lapply(iter_results[seq_len(plot_iterations)], `[[`, "transformation_curves")
        trans_curves <- trans_curves[!vapply(trans_curves, is.null, logical(1))]
        if (length(trans_curves) > 0) {
          transformation_curve_list[[paste(scenario, n_current, sep = "_")]] <- dplyr::bind_rows(trans_curves)
        }
      }
    }
  }

  metrics_df <- dplyr::bind_rows(metrics_list)
  summary_df <- summarize_simulation_metrics(metrics_df)
  action_ranges <- if (length(action_list)) dplyr::bind_rows(action_list) else NULL
  transformation_metrics <- if (length(transformation_metric_list)) dplyr::bind_rows(transformation_metric_list) else NULL
  transformation_curves <- if (length(transformation_curve_list)) dplyr::bind_rows(transformation_curve_list) else NULL

  list(
    metrics = metrics_df,
    summary = summary_df,
    action_ranges = action_ranges,
    transformation_metrics = transformation_metrics,
    transformation_curves = transformation_curves
  )
}

# -------------------- Real data utilities --------------------

get_available_real_datasets <- function() {
  c("forestfires", "bikesharing", "parkinsons")
}

get_target_features <- function(dataset_name) {
  switch(
    dataset_name,
    forestfires = c("area"),
    bikesharing = c("hum", "windspeed"),
    parkinsons = c("NHR", "HNR", "DFA", "motor_UPDRS", "total_UPDRS"),
    stop("Unknown dataset: ", dataset_name, call. = FALSE)
  )
}

load_real_dataset <- function(dataset_name, target_feature,
                              scale_y = TRUE, scale_X = TRUE) {
  dataset_name <- match.arg(dataset_name, get_available_real_datasets())

  loader <- switch(
    dataset_name,
    forestfires = function() read.csv(sebart_data_path("forestfires.csv")),
    bikesharing = function() read.csv(sebart_data_path("hour.csv")),
    parkinsons = function() read.csv(sebart_data_path("parkinsons_updrs.data"))
  )
  data <- loader()

  if (dataset_name == "parkinsons") {
    data$motor_UPDRS <- data$motor_UPDRS / 100
    data$total_UPDRS <- data$total_UPDRS / 200
  }

  if (is.null(target_feature)) {
    target_feature <- get_target_features(dataset_name)[1]
  }
  if (!target_feature %in% names(data)) {
    stop("Target feature not found in dataset: ", target_feature, call. = FALSE)
  }

  y_raw <- data[[target_feature]]
  exclude_cols <- unique(c(target_feature, "instant", "dteday", "casual",
                           "registered", "cnt"))
  X_raw <- data[, setdiff(names(data), exclude_cols), drop = FALSE]
  X_raw <- X_raw[, sapply(X_raw, is.numeric), drop = FALSE]

  idx_complete <- complete.cases(X_raw) & !is.na(y_raw)
  y_raw <- y_raw[idx_complete]
  X_raw <- X_raw[idx_complete, , drop = FALSE]

  y_mean <- mean(y_raw)
  y_sd <- sd(y_raw)
  y_min <- min(y_raw)
  y_max <- max(y_raw)

  y_proc <- if (scale_y) (y_raw - y_min) / (y_max - y_min) else y_raw
  X_proc <- if (scale_X) scale(X_raw) else as.matrix(X_raw)

  list(
    y = as.numeric(y_proc),
    X = as.matrix(X_proc),
    meta = list(
      dataset = dataset_name,
      target = target_feature,
      n_obs = length(y_proc),
      n_features = ncol(X_proc),
      y_original_range = c(y_min, y_max),
      y_mean = y_mean,
      y_sd = y_sd
    )
  )
}

create_repeated_cv_splits <- function(n, folds = 10, repeats = 10, seed = 202401) {
  set.seed(seed)
  splits <- vector("list", folds * repeats)
  counter <- 1
  for (rep in seq_len(repeats)) {
    indices <- sample.int(n)
    fold_ids <- ceiling(seq_along(indices) / ceiling(n / folds))
    fold_ids <- (fold_ids - 1) %% folds + 1
    for (fold in seq_len(folds)) {
      test_idx <- indices[fold_ids == fold]
      train_idx <- setdiff(indices, test_idx)
      splits[[counter]] <- list(
        train = sort(train_idx),
        test = sort(test_idx),
        repeat_id = rep,
        fold = fold
      )
      counter <- counter + 1
    }
  }
  splits
}

run_real_data_iteration <- function(dataset, split, models,
                                    model_config = list(ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE),
                                    alpha = 0.1,
                                    ppc_sample_size = 500,
                                    iteration_id = NULL,
                                    total_iterations = NULL,
                                    verbose_progress = FALSE) {
  y <- dataset$y
  X <- dataset$X
  train_idx <- split$train
  test_idx <- split$test

  fits <- fit_models(
    y_train = y[train_idx],
    X_train = X[train_idx, , drop = FALSE],
    X_test = X[test_idx, , drop = FALSE],
    models = models,
    config = model_config,
    progress = list(
      verbose = verbose_progress,
      iteration = iteration_id,
      total_iterations = total_iterations
    )
  )

  metrics <- do.call(rbind, lapply(names(fits), function(model_name) {
    evaluate_model(fits[[model_name]], y[test_idx], model_name, alpha)
  }))
  metrics$Repeat <- split$repeat_id
  metrics$Fold <- split$fold

  sample_records <- list()

  observed_df <- data.frame(
    Repeat = split$repeat_id,
    Fold = split$fold,
    Observation = seq_along(test_idx),
    Observed = y[test_idx],
    stringsAsFactors = FALSE
  )

  for (model_name in names(fits)) {
    draws <- fits[[model_name]]$post_ypred
    if (!is.null(draws)) {
      sampled <- sample_posterior_predictions(draws, sample_size = ppc_sample_size)
      if (!is.null(sampled)) {
        sample_records[[model_name]] <- data.frame(
          Model = model_name,
          Repeat = split$repeat_id,
          Fold = split$fold,
          Sample = sampled,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    metrics = metrics,
    posterior_samples = if (length(sample_records)) dplyr::bind_rows(sample_records) else NULL,
    observed = observed_df
  )
}

summarize_real_data_metrics <- function(metrics_df) {
  dplyr::summarise(
    dplyr::group_by(metrics_df, Dataset, Target, Model),
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    Mean_Coverage_Q_80 = mean(Coverage_Q_80, na.rm = TRUE),
    Mean_Width_Q_80 = mean(Width_Q_80, na.rm = TRUE),
    Mean_Coverage_Q_90 = mean(Coverage_Q_90, na.rm = TRUE),
    Mean_Width_Q_90 = mean(Width_Q_90, na.rm = TRUE),
    Mean_Coverage_Q_95 = mean(Coverage_Q_95, na.rm = TRUE),
    Mean_Width_Q_95 = mean(Width_Q_95, na.rm = TRUE),
    Mean_Coverage_HPD_80 = mean(Coverage_HPD_80, na.rm = TRUE),
    Mean_Width_HPD_80 = mean(Width_HPD_80, na.rm = TRUE),
    Mean_Coverage_HPD_90 = mean(Coverage_HPD_90, na.rm = TRUE),
    Mean_Width_HPD_90 = mean(Width_HPD_90, na.rm = TRUE),
    Mean_Coverage_HPD_95 = mean(Coverage_HPD_95, na.rm = TRUE),
    Mean_Width_HPD_95 = mean(Width_HPD_95, na.rm = TRUE),
    Mean_CRPS = mean(Mean_CRPS, na.rm = TRUE),
    Median_CRPS = median(Median_CRPS, na.rm = TRUE),
    Mean_Runtime = mean(Runtime, na.rm = TRUE),
    .groups = "drop"
  )
}

run_real_data_study <- function(dataset_targets,
                                models = c("sebart", "bart", "bart_bc", "sblm"),
                                folds = 10,
                                repeats = 10,
                                seed = 202402,
                                alpha = 0.1,
                                verbose = TRUE) {
  ensure_packages(c("dplyr"))
  all_metrics <- list()
  sample_list <- list()
  observed_list <- list()
  model_config <- list(ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE)

  for (dt in dataset_targets) {
    dataset_name <- dt$dataset
    target_feature <- dt$target
    if (verbose) message(sprintf("Running dataset %s (%s)", dataset_name, target_feature))

    dataset <- load_real_dataset(dataset_name, target_feature)
    splits <- create_repeated_cv_splits(dataset$meta$n_obs, folds = folds, repeats = repeats, seed = seed)

    total_splits <- length(splits)
    dt_results <- lapply(seq_along(splits), function(split_idx) {
      split <- splits[[split_idx]]
      run_real_data_iteration(
        dataset = dataset,
        split = split,
        models = models,
        model_config = model_config,
        alpha = alpha,
        ppc_sample_size = 500,
        iteration_id = split_idx,
        total_iterations = total_splits,
        verbose_progress = verbose
      )
    })

    dt_metrics <- dplyr::bind_rows(lapply(dt_results, function(res) res$metrics))
    dt_metrics$Dataset <- dataset_name
    dt_metrics$Target <- target_feature
    all_metrics[[paste(dataset_name, target_feature, sep = "_")]] <- dt_metrics

    samples <- lapply(dt_results, function(res) res$posterior_samples)
    samples <- samples[!vapply(samples, is.null, logical(1))]
    if (length(samples)) {
      sample_df <- dplyr::bind_rows(samples)
      sample_df$Dataset <- dataset_name
      sample_df$Target <- target_feature
      sample_list[[paste(dataset_name, target_feature, sep = "_")]] <- sample_df
    }

    observed_df <- dplyr::bind_rows(lapply(dt_results, function(res) res$observed))
    observed_df$Dataset <- dataset_name
    observed_df$Target <- target_feature
    observed_list[[paste(dataset_name, target_feature, sep = "_")]] <- observed_df
  }

  metrics_df <- dplyr::bind_rows(all_metrics)
  summary_df <- summarize_real_data_metrics(metrics_df)
  sample_df <- if (exists("sample_list") && length(sample_list)) dplyr::bind_rows(sample_list) else NULL
  observed_df <- if (length(observed_list)) dplyr::bind_rows(observed_list) else NULL
  list(metrics = metrics_df, summary = summary_df, posterior_samples = sample_df, observed = observed_df)
}

plot_real_data_widths <- function(metrics_df, title = NULL) {
  ensure_packages(c("ggplot2", "dplyr"))
  df <- metrics_df[!is.na(metrics_df$Width_Q_90), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  coverage <- dplyr::summarise(
    dplyr::group_by(df, Dataset, Target, Model),
    Mean_Coverage = mean(Coverage_Q_90, na.rm = TRUE),
    .groups = "drop"
  )

  df$Dataset_Target <- paste(df$Dataset, df$Target, sep = " :: ")
  coverage$Dataset_Target <- paste(coverage$Dataset, coverage$Target, sep = " :: ")
  max_width <- max(df$Width_Q_90) * 1.15

  ggplot2::ggplot(df, ggplot2::aes(x = Width_Q_90, y = Model)) +
    ggplot2::geom_boxplot(fill = "white", colour = "black", outlier.size = 1,
                          outlier.shape = 1, size = 0.5) +
    ggplot2::geom_text(
      data = coverage,
      ggplot2::aes(x = max_width, y = Model,
                   label = paste0(round(Mean_Coverage * 100), "%")),
      colour = "blue", hjust = 0, size = 3.3
    ) +
    ggplot2::facet_wrap(~Dataset_Target, ncol = 1) +
    ggplot2::labs(x = "Interval width", y = NULL,
                  title = title %||% "90% prediction interval widths") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.20)))
}

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

create_dataset_tables <- function(summary_df) {
  split(summary_df, interaction(summary_df$Dataset, summary_df$Target, drop = TRUE))
}

plot_ecdf_ppc <- function(sample_df, observed_df, dataset, target) {
  ensure_packages(c("ggplot2", "dplyr"))
  if (is.null(sample_df) || is.null(observed_df)) return(NULL)
  samples <- sample_df[sample_df$Dataset == dataset & sample_df$Target == target, , drop = FALSE]
  observed <- observed_df[observed_df$Dataset == dataset & observed_df$Target == target, , drop = FALSE]
  if (!nrow(samples) || !nrow(observed)) return(NULL)

  observed_data <- data.frame(
    Model = "Observed",
    Value = observed$Observed,
    stringsAsFactors = FALSE
  )

  posterior_data <- data.frame(
    Model = samples$Model,
    Value = samples$Sample,
    stringsAsFactors = FALSE
  )

  combined <- dplyr::bind_rows(posterior_data, observed_data)

  ggplot2::ggplot(combined, ggplot2::aes(x = Value, colour = Model)) +
    ggplot2::stat_ecdf(size = 0.9) +
    ggplot2::labs(title = sprintf("Posterior predictive ECDF: %s - %s", dataset, target),
                  x = "Value", y = "ECDF", colour = "Model") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_density_ppc <- function(sample_df, observed_df, dataset, target,
                             adjust = 1) {
  ensure_packages(c("ggplot2", "dplyr"))
  if (is.null(sample_df) || is.null(observed_df)) return(NULL)
  samples <- sample_df[sample_df$Dataset == dataset & sample_df$Target == target, , drop = FALSE]
  observed <- observed_df[observed_df$Dataset == dataset & observed_df$Target == target, , drop = FALSE]
  if (!nrow(samples) || !nrow(observed)) return(NULL)

  density_list <- lapply(split(samples, samples$Model), function(df) {
    dens <- stats::density(df$Sample, adjust = adjust)
    data.frame(Model = unique(df$Model), x = dens$x, y = dens$y, stringsAsFactors = FALSE)
  })
  density_df <- if (length(density_list)) do.call(rbind, density_list) else NULL

  obs_density <- stats::density(observed$Observed, adjust = adjust)
  obs_df <- data.frame(x = obs_density$x, y = obs_density$y)

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = density_df,
      ggplot2::aes(x = x, y = y, colour = Model),
      linewidth = 0.9
    ) +
    ggplot2::geom_line(
      data = obs_df,
      ggplot2::aes(x = x, y = y),
      colour = "black",
      linewidth = 1.1,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = sprintf("Posterior predictive densities: %s - %s", dataset, target),
      x = "Value",
      y = "Density",
      colour = "Model"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
