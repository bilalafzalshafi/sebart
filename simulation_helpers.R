#' Generate data for SBART testing scenarios using Friedman's function
#' 
#' @param n_train number of training observations
#' @param n_test number of test observations  
#' @param p number of predictors
#' @param scenario transformation scenario
#' @param heterosked logical; add heteroskedasticity
#' @param seed random seed for reproducibility
simulate_sbart_data <- function(n_train = 200, n_test = 1000, p = 10, 
                               scenario = "box_cox", heterosked = FALSE, 
                               seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  X_train <- matrix(runif(n_train * p), n_train, p)
  X_test <- matrix(runif(n_test * p), n_test, p)
  
  f_true_train <- 10 * sin(pi * X_train[,1] * X_train[,2]) + 
                  20 * (X_train[,3] - 0.5)^2 + 
                  10 * X_train[,4] + 5 * X_train[,5]
  
  f_true_test <- 10 * sin(pi * X_test[,1] * X_test[,2]) + 
                 20 * (X_test[,3] - 0.5)^2 + 
                 10 * X_test[,4] + 5 * X_test[,5]
  
  if(heterosked) {
    z_train <- f_true_train * (1 + rnorm(n_train, sd = 0.5))
    z_test <- f_true_test * (1 + rnorm(n_test, sd = 0.5))
  } else {
    z_train <- f_true_train + rnorm(n_train, sd = 0.5)
    z_test <- f_true_test + rnorm(n_test, sd = 0.5)
  }
  
  # Apply g_inv to latent z to get observed y (like simulate_tlm)
  transform_info <- get_transformation(scenario)
  y_train <- transform_info$g_inv(z_train)
  y_test <- transform_info$g_inv(z_test)
  
  # True transformation g for diagnostics
  y_unique <- sort(unique(y_train))
  g_true <- if(is.null(transform_info$g)) NULL else transform_info$g(y_unique)
  
  list(
    X_train = X_train, y_train = y_train,
    X_test = X_test, y_test = y_test,
    f_true_train = f_true_train, f_true_test = f_true_test,
    z_train = z_train, z_test = z_test,
    scenario = scenario,
    g_inverse = transform_info$g_inv,
    g_true = g_true,
    y_unique = y_unique,
    description = transform_info$description
  )
}

get_transformation <- function(scenario) {
  transformations <- list(

    "step" = list(
      g = function(y) {y},
      g_inv = function(z) {
        # Create monotonic step function
        g_steps = rexp(n = 10)
        approxfun(seq(-3, 3, length.out=10),
                  cumsum(g_steps), rule = 2)((z - mean(z)) / sd(z))
      },
      description = "Monotonic step transformation"
    ),

    "sigmoid" = list(
      g = function(y) -2 * log(20/pmax(y, 1e-10) - 1),
      g_inv = function(z) 20 / (1 + exp(-z/2)),
      description = "Inverse sigmoid transformation, bounded [0,20]"
    ),

    "beta" = list(
      g = function(y) qnorm(pbeta(pmax(pmin(y, 1-1e-10), 1e-10), 0.1, 0.5)),
      g_inv = function(z) {
        z_std <- (z - mean(z)) / sd(z)
        qbeta(pnorm(z_std), shape1 = 0.1, shape2 = 0.5)
      },
      description = "Beta(0.1,0.5) bounded [0,1], many values near zero"
    ),

    "arctangent" = list(
      g = function(y) 5 * tan(y/100),
      g_inv = function(z) 100 * atan(z/5),
      description = "Arctangent transformation, smooth bounded"
    ),

    "box_cox" = list(
      g = function(y) g_bc(y, lambda = 0.5),
      g_inv = function(z) SebR:::g_inv_bc(z, lambda = 0.5),
      description = "Inverse Box-Cox square root transformation"
    ),

    "almost_linear" = list(
      g = function(y) (-1 + sqrt(1 + 0.4*y)) / 0.2,
      g_inv = function(z) z + 0.1 * z^2,
      description = "Nearly linear with weak quadratic component"
    ),

    "polynomial" = list(
      g = function(y) y^(1/2.5),
      g_inv = function(z) {
        z_positive <- z - min(z) + 5
        z_positive^2.5
      },
      description = "Polynomial x^2.5 transformation"
    )
  )

  if (!scenario %in% names(transformations)) {
    stop("Unknown scenario. Available: ", paste(names(transformations), collapse = ", "))
  }

  transformations[[scenario]]
}

get_available_scenarios <- function() {
  c("step", "sigmoid", "beta", "arctangent", "box_cox", "almost_linear", "polynomial")
}