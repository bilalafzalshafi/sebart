#' Semiparametric Bayesian BART with Location-Scale Robustness
#'
#' MCMC sampling for semiparametric Bayesian regression with
#' 1) an unknown (nonparametric) transformation and 2) a BART model
#' for the regression function. The transformation is modeled as unknown
#' and learned jointly with the BART model using a blocked Gibbs sampler
#' with proper location-scale parameter expansion for robustness.
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
#' \item \code{post_mu}: \code{nsave} posterior samples of the location parameter
#' \item \code{post_sigma}: \code{nsave} posterior samples of the scale parameter
#' \item \code{model}: the model fit (here, \code{sbart})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed BART model using blocked Gibbs sampling with proper
#' location-scale parameter expansion (LS-PX) for robustness.
#' The transformation is modeled as unknown and learned jointly
#' with the BART model (unless \code{approx_g} = TRUE, which then uses
#' a point approximation). This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' By default, \code{fixedX} is set to \code{FALSE} for smaller datasets (\code{n < 500})
#' and \code{TRUE} for larger datasets.
#'
#' The approach uses a pilot BART run to approximate the posterior distribution
#' required for transformation inference. The main MCMC then alternates
#' between sampling the transformation (given current BART fit), sampling
#' location-scale parameters (μ,σ) for robustness, and updating the BART model.
#'
#' @note The location-scale parameters (μ,σ) are properly sampled from their
#' posterior distribution and the reported transformation is adjusted to ensure
#' identifiability following the LS-PX framework.
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
  
  # Hyperparameters for Gamma(a_sigma, b_sigma) prior on error precision
  a_sigma = b_sigma = 0.001
  
  #----------------------------------------------------------------------------
  # Initialize the transformation:

  # Define the CDF of y:
  Fy = function(t) n/(n+1)*ecdf(y)(t)

  # Evaluate at the unique y-values:
  y0 = sort(unique(y))
  Fy_eval = Fy(y0)

  # Initial transformation: quantile-quantile transformation
  z = qnorm(Fy(y))

  #----------------------------------------------------------------------------
  # Pilot BART run to approximate posterior
  if(verbose) cat('Running pilot BART to approximate posterior...\n')
  
  # Construct data frame with matrix columns
  X_df <- as.data.frame(X)
  colnames(X_df) <- paste0("X", 1:ncol(X))
  pilot_data <- cbind(data.frame(z = z), X_df)
  
  # Set up pilot BART sampler
  pilot_sampler = dbarts::dbarts(
    formula = z ~ .,
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
  pilot_sigma = mean(pilot_samples$sigma)
  
  if(verbose) cat('Pilot run complete. Initializing main sampler...\n')

  #----------------------------------------------------------------------------
  # Grid of values for the CDF of z based on pilot estimates:
  z_grid = sort(unique(
    sapply(range(pilot_fits), function(xtemp){
      qnorm(seq(0.01, 0.99, length.out = ngrid),
            mean = xtemp,
            sd = pilot_sigma)
    })
  ))

  # Store F_{Z | X = x_i}(t) for all x_i & all t in z_grid
  # This integrates out the BART function from its posterior:
  zrep = rep(z_grid, times = pilot_ndraws)
  Fzx_eval = sapply(1:n, function(i){
    rowMeans(matrix(pnorm(zrep,
                          mean = rep(pilot_fits[i,], each = ngrid),
                          sd = pilot_sigma), nrow = ngrid))
  })

  # CDF of z using pilot estimates (initial approximation):
  Fz_eval = rowMeans(Fzx_eval)  # average across observations for initial grid

  # Check: update the grid if needed
  zcon = contract_grid(z = z_grid,
                       Fz = Fz_eval,
                       lower = 0.001, upper = 0.999)
  z_grid = zcon$z; Fz_eval = zcon$Fz; ngrid = length(z_grid)
  
  # Recompute Fzx_eval with updated grid
  zrep = rep(z_grid, times = pilot_ndraws)
  Fzx_eval = sapply(1:n, function(i){
    rowMeans(matrix(pnorm(zrep,
                          mean = rep(pilot_fits[i,], each = ngrid),
                          sd = pilot_sigma), nrow = ngrid))
  })

  # Recompute Fz_eval with updated grid
  Fz_eval = rowMeans(Fzx_eval)

  # Compute initial transformation:
  g = g_fun(y = y0,
            Fy_eval = Fy_eval,
            z = z_grid,
            Fz_eval = Fz_eval)

  # Apply initial transformation
  z = g(y)

  # Define the grid for approximations using equally-spaced + quantile points:
  y_grid = sort(unique(c(
    seq(min(y), max(y), length.out = ngrid/2),
    quantile(y0, seq(0, 1, length.out = ngrid/2)))))

  # Inverse transformation function:
  g_inv = g_inv_approx(g = g, t_grid = y_grid)

  #----------------------------------------------------------------------------
  # Initialize location-scale parameters following LS-PX framework
  mu = mean(z)  # location parameter
  sigma_epsilon = sd(z)  # scale parameter
  
  # Standardize z for BART (following LS-PX: z_standardized = (z - mu)/sigma)
  z_standardized = (z - mu)/sigma_epsilon
  
  # Set up main BART sampler with standardized data
  X_df <- as.data.frame(X)
  colnames(X_df) <- paste0("X", 1:ncol(X))
  main_data <- cbind(data.frame(z_std = z_standardized), X_df)
  
  # # Set up main sampler
  # main_sampler = dbarts::dbarts(
  #   formula = z_std ~ .,
  #   data = main_data,
  #   control = dbarts::dbartsControl(
  #     verbose = FALSE,
  #     keepTrainingFits = TRUE,
  #     keepTrees = FALSE,
  #     n.trees = ntree,
  #     n.chains = 1L,
  #     n.threads = n.threads,
  #     rngSeed = seed,
  #     updateState = FALSE
  #   )
  # )

  # Calculate properties of double-scaled data for prior calibration
  z_std_var = var(z_standardized)
  z_std_range = diff(range(z_standardized))

  # Estimate what BART's internal rescaling will produce
  # BART rescales to [-0.5, 0.5], so final range will be 1.0
  # Final variance will be approximately z_std_var * (1.0/z_std_range)^2
  final_expected_var = z_std_var * (1.0/z_std_range)^2

  # Adjust k: default k=2 assumes certain scale properties
  # If our double-scaled data has different variance, adjust proportionally
  k_adjusted = 2.0 * sqrt(final_expected_var / 0.25)  # 0.25 is var for uniform[-0.5,0.5]

  # Set up main sampler with calibrated priors
  main_sampler = dbarts::dbarts(
    formula = z_std ~ .,
    data = main_data,
    sigma = sqrt(final_expected_var),
    control = dbarts::dbartsControl(
      verbose = FALSE,
      keepTrainingFits = TRUE,
      keepTrees = FALSE,
      n.trees = ntree,
      n.chains = 1L,
      n.threads = n.threads,
      rngSeed = seed,
      updateState = FALSE,
    ),
    node.prior = normal(k = k_adjusted),
    resid.prior = chisq(df = 3, quant = 0.90)
  )

  # Initialize MCMC sampler using object method
  invisible(main_sampler$run(numBurnIn = 1, numSamples = 1))

  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))
  post_mu = rep(NA, nsave)
  post_sigma = rep(NA, nsave)

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

        # CDF of Z evaluated on the grid using pilot posterior:
        Fz_eval = Fzx_eval%*%weights_x
      } else {
        # Fixed X case - use average across observations
        Fz_eval = rowMeans(Fzx_eval)
      }

      # Compute the transformation:
      g = g_fun(y = y0, Fy_eval = Fy_eval,
                z = z_grid, Fz_eval = Fz_eval)

      # Update z:
      z = g(y)

      # Update the inverse transformation function:
      g_inv = g_inv_approx(g = g, t_grid = y_grid)
    }

    #----------------------------------------------------------------------------
    # Block 2: sample location-scale parameters (μ,σ) following LS-PX framework
    
    # Get current BART fit for residual calculation
    current_sample = main_sampler$run(numBurnIn = 0, numSamples = 1)
    current_fit_std = as.vector(current_sample$train) # standardized scale
    current_fit = mu + sigma_epsilon * current_fit_std # original scale
    
    # Residuals on original scale
    residuals = z - current_fit
    
    # Sample sigma_epsilon (scale parameter) using conjugate gamma prior
    SSR = sum(residuals^2)
    sigma_epsilon = 1/sqrt(rgamma(n = 1,
                                  shape = a_sigma + n/2,
                                  rate = b_sigma + SSR/2))
    
    # Sample mu (location parameter) using conjugate normal prior
    # Prior: mu ~ N(0, large_var), which is essentially flat
    mu_var = 1000  # diffuse prior variance
    mu_post_var = 1/(n/sigma_epsilon^2 + 1/mu_var)
    mu_post_mean = mu_post_var * sum(z - current_fit + mu)/sigma_epsilon^2
    mu = rnorm(1, mean = mu_post_mean, sd = sqrt(mu_post_var))
    
    # Update standardized z for BART
    z_standardized = (z - mu)/sigma_epsilon
    
    # Update BART sampler with new standardized response
    main_sampler$setResponse(z_standardized)

    #----------------------------------------------------------------------------
    # Block 3: sample BART model parameters and functions
    
    # Single MCMC step for BART using object method
    bart_sample = main_sampler$run(numBurnIn = 0, numSamples = 1)

    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Predictive samples:
      # Get BART predictions at test points using object method
      X_test_df <- as.data.frame(X_test)
      colnames(X_test_df) <- paste0("X", 1:ncol(X_test))
      bart_pred_std = main_sampler$predict(X_test_df)
      
      # Transform back to original scale: z_tilde = mu + sigma * f_BART(x) + sigma * epsilon
      current_sigma_std = bart_sample$sigma[1]
      bart_pred = mu + sigma_epsilon * as.vector(bart_pred_std)
      ztilde = bart_pred + sigma_epsilon * current_sigma_std * rnorm(n = n_test)
      post_ypred[nsi - nburn,] = g_inv(ztilde)

      # Posterior samples of location-scale parameters:
      post_mu[nsi - nburn] = mu
      post_sigma[nsi - nburn] = sigma_epsilon

      # Posterior samples of the transformation (adjusted for location-scale following LS-PX):
      # Report g^(μ,σ) = (g(y0) - μ)/σ for identifiability
      post_g[nsi - nburn,] = (g(y0) - mu)/sigma_epsilon
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
    post_mu = post_mu,
    post_sigma = post_sigma,
    model = 'sbart', 
    y = y, X = X, X_test = X_test, 
    fixedX = fixedX, approx_g = approx_g,
    pilot_ndraws = pilot_ndraws, pilot_nburn = pilot_nburn,
    ntree = ntree,
    n.threads = n.threads, seed = seed))
}