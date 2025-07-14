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
  
  # Generate design matrices
  X_train <- matrix(runif(n_train * p), n_train, p)
  X_test <- matrix(runif(n_test * p), n_test, p)
  
  # Friedman's nonlinear function
  f_true_train <- 10 * sin(pi * X_train[,1] * X_train[,2]) + 
                  20 * (X_train[,3] - 0.5)^2 + 
                  10 * X_train[,4] + 5 * X_train[,5]
  
  f_true_test <- 10 * sin(pi * X_test[,1] * X_test[,2]) + 
                 20 * (X_test[,3] - 0.5)^2 + 
                 10 * X_test[,4] + 5 * X_test[,5]
  
  # Add noise
  if(heterosked) {
    z_train <- f_true_train * (1 + rnorm(n_train))
    z_test <- f_true_test * (1 + rnorm(n_test))
  } else {
    z_train <- f_true_train + rnorm(n_train, sd = 0.5)
    z_test <- f_true_test + rnorm(n_test, sd = 0.5)
  }
  
  # Standardize for some transformations (matching diagnostics logic)
  if (scenario %in% c("beta", "sigmoid")) {
    z_train <- (z_train - mean(z_train)) / sd(z_train)
    z_test <- (z_test - mean(c(z_train, z_test))) / sd(c(z_train, z_test))
  }
  
  
  # Apply transformation
  transform_info <- get_transformation(scenario)
  y_train <- transform_info$transform(z_train)
  y_test <- transform_info$transform(z_test)
  
  # Compute true transformation values for diagnostics
  y_unique <- sort(unique(y_train))
  g_true <- NULL
  if (!is.null(transform_info$inverse)) {
    g_true <- transform_info$inverse(y_unique)
  }
  
  list(
    X_train = X_train, y_train = y_train,
    X_test = X_test, y_test = y_test,
    f_true_train = f_true_train, f_true_test = f_true_test,
    z_train = z_train, z_test = z_test,
    scenario = scenario,
    g_inverse = transform_info$inverse,
    g_true = g_true,
    y_unique = y_unique,
    description = transform_info$description
  )
}

#' Get transformation function and metadata
get_transformation <- function(scenario) {
  transformations <- list(
    "box_cox" = list(
      transform = function(z) {
        z_shifted <- z - min(z) + 1  # Ensure positive
        z_shifted^2  # Square transformation
      },
      inverse = function(y) sqrt(pmax(y, 0)),
      description = "Box-Cox lambda=0.5 (square root)"
    ),
    
    "step" = list(
      transform = function(z) {
        q30 <- quantile(z, 0.3)
        q70 <- quantile(z, 0.7)
        ifelse(z < q30, z^2, 
               ifelse(z < q70, 2*z + 5, 
                      exp((z - q70)*0.1) + 2*q70 + 5))
      },
      inverse = NULL,  # Complex piecewise inverse
      description = "Step function, positive support"
    ),
    
    "sigmoid" = list(
      transform = function(z) 20 / (1 + exp(-z/2)),
      inverse = function(y) -2 * log(20/y - 1),
      description = "Sigmoid, bounded [0,20]"
    ),
    
    "linear" = list(
      transform = function(z) z + 0.1 * z^2,
      inverse = function(y) (-1 + sqrt(1 + 0.4*y)) / 0.2,  # Quadratic formula
      description = "Nearly linear with weak quadratic"
    ),
    
    "identity" = list(
      transform = function(z) z,
      inverse = function(y) y,
      description = "Identity (no transformation)"
    ),
    
    "beta" = list(
      transform = function(z) {
        z_std <- (z - mean(z)) / sd(z)
        qbeta(pnorm(z_std), shape1 = 0.1, shape2 = 0.5)
      },
      inverse = function(y) qnorm(pbeta(y, 0.1, 0.5)),
      description = "Beta(0.1,0.5) bounded [0,1], many values near zero"
    )
  )
  
  if (!scenario %in% names(transformations)) {
    stop("Unknown scenario. Available: ", paste(names(transformations), collapse = ", "))
  }
  
  transformations[[scenario]]
}

#' Get available transformation scenarios
get_available_scenarios <- function() {
  c("box_cox", "step", "sigmoid", "linear", "identity", "beta")
}