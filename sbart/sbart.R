#' Semiparametric Bayesian BART
#'
#' MCMC sampling for semiparametric Bayesian regression with
#' 1) an unknown (nonparametric) transformation and 2) a BART model
#' for the regression function. The transformation is modeled as unknown
#' and learned jointly with the BART model using a blocked Gibbs sampler.
#' A pilot BART run is used to approximate the posterior distribution
#' needed for transformation inference.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param approx_g logical; if TRUE, apply large-sample
#' approximation for the transformation
#' @param pilot_ndraws number of MCMC draws for pilot BART run to approximate posterior
#' @param pilot_nburn number of burn-in iterations for pilot BART run  
#' @param ntree number of trees in BART ensemble
#' @param nsave number of MCMC simulations to save
#' @param nburn number of MCMC iterations to discard
#' @param ngrid number of grid points for inverse approximations
#' @param verbose logical; if TRUE, print time remaining
#' @param n.threads integer; number of threads to use (default 1)
#' @param seed integer; random seed for reproducibility (default NA)
#' @return a list with the following elements:
#' \itemize{
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{X_test}
#' \item \code{post_ypred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{post_sigma}: \code{nsave} posterior samples of the error standard deviation
#' \item \code{model}: the model fit (here, \code{sbart})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed BART model using blocked Gibbs sampling.
#' The transformation is modeled as unknown and learned jointly
#' with the BART model (unless \code{approx_g} = TRUE, which then uses
#' a point approximation). This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' By default, \code{fixedX} is set to \code{FALSE} for smaller datasets (\code{n < 500})
#' and \code{TRUE} for larger datasets.
#'
#' The approach uses a pilot BART run to approximate the posterior distribution
#' required for transformation inference. This pilot run provides estimates of
#' E[f_BART(x)|data] and Var[f_BART(x)|data] which are used to approximate
#' the CDF F_Z needed for transformation sampling. The main MCMC then alternates
#' between sampling the transformation (given current BART fit) and updating
#' the BART model (given current transformation).
#'
#' @note The location (intercept) and scale are not identified under transformation
#' models, so the BART model includes internal location-scale adjustments but
#' the function outputs only refer to the identifiable transformed scale.
#'
#' @examples
#' \donttest{
#' # Simulate data from a transformed nonlinear model:
#' n <- 200
#' p <- 5
#' X <- matrix(rnorm(n*p), n, p)
#' f_true <- sin(2*X[,1]) + X[,2]^2 - X[,3]*X[,4] + rnorm(n, sd=0.1)
#' y <- exp(f_true + 0.5*rnorm(n))  # log-normal with nonlinear mean
#' 
#' # Test data
#' X_test <- matrix(rnorm(50*p), 50, p)
#' 
#' # Fit the semiparametric Bayesian BART model:
#' fit <- sbart(y = y, X = X, X_test = X_test)
#' names(fit) # what is returned
#' 
#' # Evaluate posterior predictive intervals on testing data:
#' pi_y <- t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
#' plot(1:nrow(pi_y), pi_y[,1], type='n', ylim=range(pi_y),
#'      xlab='Test observation', ylab='y', main='Prediction intervals')
#' segments(1:nrow(pi_y), pi_y[,1], 1:nrow(pi_y), pi_y[,2])
#' points(1:nrow(pi_y), fitted(fit), pch=16)
#' 
#' # Summarize the transformation:
#' y0 <- sort(unique(y)) # posterior draws of g are evaluated at unique y values
#' plot(y0, fit$post_g[1,], type='n', ylim=range(fit$post_g),
#'      xlab='y', ylab='g(y)', main="Posterior draws of the transformation")
#' for(s in sample(nrow(fit$post_g), 50)) lines(y0, fit$post_g[s,], col='gray')
#' lines(y0, colMeans(fit$post_g), lwd=3) # posterior mean
#' }
#' @export
sbart = function(y, X, X_test = X,
                 fixedX = (length(y) >= 500),
                 approx_g = FALSE,
                 pilot_ndraws = 500,
                 pilot_nburn = 100,
                 ntree = 200,
                 nsave = 1000,
                 nburn = 1000,
                 ngrid = 100,
                 verbose = TRUE,
                 n.threads = 1,
                 seed = NA_integer_){

  # Library required here:
  if (!requireNamespace("dbarts", quietly = TRUE)) {
    stop(
      "Package \"dbarts\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Initial checks:
  if(!is.matrix(X)) stop("X must be a matrix (rows = observations, columns = variables)")
  if(!is.matrix(X_test)) stop("X_test must be a matrix (rows = observations, columns = variables)")
  if(ncol(X) != ncol(X_test)) stop('X_test and X must have the same number of columns (variables)')
  if(nrow(X) != length(y)) stop('the length of y must equal nrow(X)')
  if(any(is.na(y)) || any(is.na(X))) stop('y and X cannot contain missing values')
  if(any(!is.finite(y)) || any(!is.finite(X))) stop('y and X must be finite')

  # Data dimensions:
  n = length(y) # number of observations
  p = ncol(X) # number of variables
  n_test = nrow(X_test) # number of testing data points
  
  # Parameter validation
  if(pilot_ndraws < 50) stop('pilot_ndraws must be at least 50')
  if(ntree < 10) stop('ntree must be at least 10')
  
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)

  # Evaluate at the unique y-values:
  y0 = sort(unique(y))
  Fy_eval = Fy(y0)

  # Initial transformation: quantile-quantile transformation
  z = qnorm(Fy(y))
  z_scaled = (z - mean(z))/sd(z) # standardize for BART

  #----------------------------------------------------------------------------
  # Pilot BART run to approximate posterior
  if(verbose) cat('Running pilot BART to approximate posterior...\n')
  
  # Construct data frame with matrix columns
  X_df <- as.data.frame(X)
  colnames(X_df) <- paste0("X", 1:ncol(X))
  pilot_data <- cbind(data.frame(z_scaled = z_scaled), X_df)
  
  # Set up pilot BART sampler
  pilot_sampler = dbarts::dbarts(
    formula = z_scaled ~ .,
    data = pilot_data,
    control = dbarts::dbartsControl(
      verbose = FALSE,
      keepTrainingFits = TRUE,
      n.trees = ntree,
      n.chains = 1L,
      n.threads = 1L,
      rngSeed = seed
    )
  )

  # Run pilot MCMC using object method
  pilot_samples = pilot_sampler$run(numBurnIn = pilot_nburn, numSamples = pilot_ndraws)

  # Extract posterior summaries from pilot run
  pilot_fits = pilot_samples$train # n x pilot_ndraws matrix
  mu_z_pilot = rowMeans(pilot_fits) # posterior mean E[f(xi) | data]
  var_z_pilot = apply(pilot_fits, 1, var) # posterior variance Var[f(xi) | data]
  
  # Add error variance estimate from pilot
  pilot_sigma = mean(pilot_samples$sigma)
  sigma_z_pilot = sqrt(var_z_pilot + pilot_sigma^2)

  if(verbose) cat('Pilot run complete. Initializing main sampler...\n')

  #----------------------------------------------------------------------------
  # Grid of values for the CDF of z based on pilot estimates:
  z_grid = sort(unique(
    sapply(range(mu_z_pilot), function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = ngrid),
            mean = xtemp,
            sd = pilot_sigma)
    })
  ))

  # CDF of z using pilot estimates:
  Fz_eval = Fz_fun(z = z_grid,
                   weights = rep(1/n, n),
                   mean_vec = mu_z_pilot,
                   sd_vec = sigma_z_pilot)

  # Check: update the grid if needed
  zcon = contract_grid(z = z_grid,
                       Fz = Fz_eval,
                       lower = 0.001, upper = 0.999)
  z_grid = zcon$z; Fz_eval = zcon$Fz; ngrid = length(z_grid)

  # Compute initial transformation:
  g = g_fun(y = y0,
            Fy_eval = Fy_eval,
            z = z_grid,
            Fz_eval = Fz_eval)

  # Apply initial transformation
  z = g(y)
  z_scaled = (z - mean(z))/sd(z)

  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)

  #----------------------------------------------------------------------------
  # Set up main BART sampler with updated data
  # Properly construct data frame with matrix columns
  X_df <- as.data.frame(X)
  colnames(X_df) <- paste0("X", 1:ncol(X))
  main_data <- cbind(data.frame(z_scaled = z_scaled), X_df)
  
  # Set up main sampler (use defaults for priors)
  main_sampler = dbarts::dbarts(
    formula = z_scaled ~ .,
    data = main_data,
    control = dbarts::dbartsControl(
      verbose = FALSE,
      keepTrainingFits = TRUE,
      keepTrees = FALSE,
      n.trees = ntree,
      n.chains = 1L,
      n.threads = n.threads,
      rngSeed = seed,
      updateState = FALSE
    )
  )

  # Initialize MCMC sampler using object method
  invisible(main_sampler$run(numBurnIn = 1, numSamples = 1))

  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))
  post_sigma = rep(NA, nsave)

  # MCMC storage for updating transformation
  current_mu_z = mu_z_pilot
  current_sigma_z = sigma_z_pilot
  z_scale_mean = mean(z)
  z_scale_sd = sd(z)

  # Run the main MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:(nburn + nsave)){

    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(!approx_g){

      # Bayesian bootstrap for the CDFs
      # BB CDF of y:
      Fy_eval = bb(y)(y0) # applied to y, evaluated at y0

      # BB CDF of z (only need if X is random)
      if(!fixedX){

        # Dirichlet(1) weights for x:
        weights_x = rgamma(n = n, shape = 1)
        weights_x = weights_x/sum(weights_x)

        # CDF of Z evaluated on the grid using current BART estimates:
        Fz_eval = Fz_fun(z = z_grid,
                         weights = weights_x,
                         mean_vec = current_mu_z,
                         sd_vec = current_sigma_z)
      } else {
        # Fixed X case - use uniform weights
        Fz_eval = Fz_fun(z = z_grid,
                         weights = rep(1/n, n),
                         mean_vec = current_mu_z,
                         sd_vec = current_sigma_z)
      }

      # Compute the transformation:
      g = g_fun(y = y0, Fy_eval = Fy_eval,
                z = z_grid, Fz_eval = Fz_eval)

      # Update z and rescale:
      z = g(y)
      z_scale_mean = mean(z)
      z_scale_sd = sd(z)
      z_scaled = (z - z_scale_mean)/z_scale_sd

      # Update BART sampler with new response using object method:
      main_sampler$setResponse(z_scaled)

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = y_grid)
    }

    #----------------------------------------------------------------------------
    # Block 2: sample BART model parameters and functions
    
    # Single MCMC step for BART using object method
    bart_sample = main_sampler$run(numBurnIn = 0, numSamples = 1)
    
    # Extract current fit and sigma
    current_fit = as.vector(bart_sample$train) # current f(xi) estimates (scaled)
    current_sigma_scaled = bart_sample$sigma[1] # error SD (scaled)
    
    # Transform back to original scale for transformation updates
    current_mu_z = current_fit * z_scale_sd + z_scale_mean
    current_sigma_z = rep(current_sigma_scaled * z_scale_sd, n)

    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Predictive samples:
      # Get BART predictions at test points using object method
      X_test_df <- as.data.frame(X_test)
      colnames(X_test_df) <- paste0("X", 1:ncol(X_test))
      bart_pred_scaled = main_sampler$predict(X_test_df)
      bart_pred = as.vector(bart_pred_scaled) * z_scale_sd + z_scale_mean
      
      # Add noise and transform back to y scale
      ztilde = bart_pred + current_sigma_scaled * z_scale_sd * rnorm(n = n_test)
      post_ypred[nsi - nburn,] = g_inv(ztilde)

      # Posterior samples of the transformation (adjusted for location/scale):
      post_g[nsi - nburn,] = (g(y0) - z_scale_mean)/z_scale_sd

      # Posterior samples of sigma (on original scale):
      post_sigma[nsi - nburn] = current_sigma_scaled * z_scale_sd
    }

    if(verbose) computeTimeRemaining(nsi, timer0, nsave + nburn)
  }

  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer0
    if(tfinal > 60){
      cat(paste('Total time:', round(tfinal/60,1), 'minutes\n'))
    } else cat(paste('Total time:', round(tfinal), 'seconds\n'))
  }

  return(list(
    fitted.values = colMeans(post_ypred),
    post_ypred = post_ypred,
    post_g = post_g,
    post_sigma = post_sigma,
    model = 'sbart', 
    y = y, X = X, X_test = X_test, 
    fixedX = fixedX, approx_g = approx_g,
    pilot_ndraws = pilot_ndraws, pilot_nburn = pilot_nburn,
    ntree = ntree,
    n.threads = n.threads, seed = seed))
}