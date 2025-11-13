#' Semiparametric BART

#' MCMC sampling for semiparametric Bayesian regression with
#' 1) an unknown (nonparametric) transformation and 2) a BART model
#' for the regression function.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param fixedX logical; if TRUE, treat the design as fixed (non-random) when sampling
#' the transformation; otherwise treat covariates as random with an unknown distribution
#' @param pilot_method character; method for pilot approximation. Either "bart" (default)
#' or "rf" (Random Forest)
#' @param pilot_ndraws number of MCMC draws for pilot BART run to approximate posterior
#' (ignored if pilot_method = "rf")
#' @param pilot_nburn number of burn-in iterations for pilot BART run
#' (ignored if pilot_method = "rf")
#' @param rf_ntree number of trees in Random Forest (used if pilot_method = "rf")
#' @param rf_nodesize minimum size of terminal nodes in Random Forest (used if pilot_method = "rf")
#' @param rf_mtry number of variables randomly sampled as candidates at each split in Random Forest
#' (used if pilot_method = "rf"; default is max(floor(p/3), 1))
#' @param ntree number of trees in BART ensemble
#' @param nsave number of MCMC simulations to save
#' @param nburn number of MCMC iterations to discard
#' @param ngrid number of grid points for inverse approximations
#' @param verbose logical; if TRUE, print time remaining
#' @param n.threads integer; number of threads to use (default 1)
#' @param seed integer; random seed for reproducibility (default NA)
#' @param var_select logical; if TRUE, perform variable selection (default FALSE)
#' @param var_select_threshold numeric; threshold for variable selection (default 0.01)
#' @return a list with the following elements:
#' \itemize{
#' \item \code{fitted.values} the posterior predictive mean at the test points \code{X_test}
#' \item \code{post_ypred}: \code{nsave x n_test} samples
#' from the posterior predictive distribution at test points \code{X_test}
#' \item \code{post_g}: \code{nsave} posterior samples of the transformation
#' evaluated at the unique \code{y} values
#' \item \code{model}: the model fit (here, \code{sebart})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides fully Bayesian inference for a
#' transformed BART model using blocked Gibbs sampling.
#' The transformation is modeled as unknown and learned jointly
#' with the BART model. This model applies for real-valued data, positive data, and
#' compactly-supported data (the support is automatically deduced from the observed \code{y} values).
#' By default, \code{fixedX} is set to \code{FALSE} for smaller datasets (\code{n < 500})
#' and \code{TRUE} for larger datasets.
#'
#' The approach uses either a pilot BART run or Random Forest to approximate the posterior distribution
#' required for transformation inference. The main MCMC then alternates
#' between sampling the transformation (given current BART fit) on the identified
#' scale and updating the BART model, which absorbs affine adjustments through its
#' internal intercept and residual variance.
#' @export
sebart = function(y, X, X_test = X,
                 fixedX = (length(y) >= 500),
                 pilot_method = c("bart", "rf"),
                 pilot_ndraws = 500,
                 pilot_nburn = 100,
                 rf_ntree = 500,
                 rf_nodesize = 5,
                 rf_mtry = NULL,
                 ntree = 200,
                 nsave = 1000,
                 nburn = 1000,
                 ngrid = 100,
                  verbose = TRUE,
                  n.threads = 1,
                  seed = NA_integer_,
                  var_select = FALSE,
                  var_select_threshold = 0.01){

  # Match pilot method argument
  pilot_method = match.arg(pilot_method)

  # Library required here:
  if (!requireNamespace("dbarts", quietly = TRUE)) {
    stop(
      "Package \"dbarts\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  if (pilot_method == "rf" && !requireNamespace("randomForest", quietly = TRUE)) {
    stop(
      "Package \"randomForest\" must be installed to use pilot_method = 'rf'.",
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
  if(pilot_method == "bart" && pilot_ndraws < 50) stop('pilot_ndraws must be at least 50')
  if(pilot_method == "rf" && is.null(rf_mtry)) rf_mtry = max(floor(p/3), 1)
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

  #----------------------------------------------------------------------------
  # Pilot approximation to approximate posterior
  
  if(pilot_method == "bart") {
    
    if(verbose) cat('Running pilot BART to approximate posterior...\n')
    
    # Construct data frame with matrix columns
    X_df <- as.data.frame(X)
    colnames(X_df) <- paste0("X", 1:ncol(X))
    x0_df <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(X_df)))
    colnames(x0_df) <- colnames(X_df)
    pilot_data <- cbind(data.frame(z = z), X_df)
    
    # Set up pilot BART sampler
    pilot_sampler = dbarts::dbarts(
      formula = z ~ .,
      data = pilot_data,
      test = x0_df,
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
    f0_draws <- as.numeric(pilot_samples$test)  # draws of f_theta(0)
    
    if(verbose) cat('Pilot run complete. Initializing main sampler...\n')

    #--------------------------------------------------------------------------
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
    zrep_std = (zrep - rep(f0_draws, each = ngrid)) / pilot_sigma
    Fzx_eval = sapply(1:n, function(i){
      rowMeans(matrix(pnorm(zrep_std,
                            mean = rep((pilot_fits[i,] - f0_draws)/pilot_sigma,
                                       each = ngrid),
                            sd = 1), nrow = ngrid))
    })

    # CDF of z using pilot estimates (initial approximation):
    Fz_eval = rowMeans(Fzx_eval)  # average across observations for initial grid

    # Check: update the grid if needed
    zcon = SeBR:::contract_grid(z = z_grid,
                         Fz = Fz_eval,
                         lower = 0.001, upper = 0.999)
    z_grid = zcon$z; Fz_eval = zcon$Fz; ngrid = length(z_grid)
    
    # Recompute Fzx_eval with updated grid
    zrep = rep(z_grid, times = pilot_ndraws)
    zrep_std = (zrep - rep(f0_draws, each = ngrid)) / pilot_sigma
    Fzx_eval = sapply(1:n, function(i){
      rowMeans(matrix(pnorm(zrep_std,
                            mean = rep((pilot_fits[i,] - f0_draws)/pilot_sigma,
                                       each = ngrid),
                            sd = 1), nrow = ngrid))
    })

    # Recompute Fz_eval with updated grid
    Fz_eval = rowMeans(Fzx_eval)
    
  } else if(pilot_method == "rf") {
    
    if(verbose) cat("Fitting Random Forest for initial approximation...\n")
    
    # Load randomForest library
    library(randomForest)
    
    # Fit Random Forest with variance estimation
    rf_fit = randomForest(x = X, y = z, 
                          ntree = rf_ntree,
                          nodesize = rf_nodesize,
                          mtry = rf_mtry,
                          keep.inbag = TRUE,
                          importance = FALSE)
    
    # Get predictions and individual tree predictions
    rf_pred = predict(rf_fit, X, predict.all = TRUE)
    rf_mean = rf_pred$aggregate
    rf_trees = rf_pred$individual
    
    # Compute ensemble variance (uncertainty from trees)
    rf_var_ensemble = apply(rf_trees, 1, var)
    
    # Estimate residual variance using out-of-bag predictions
    oob_residuals = z - rf_fit$predicted
    sigma_epsilon_hat = sd(oob_residuals)
    
    # Total predictive variance: ensemble + residual
    rf_var_total = rf_var_ensemble + sigma_epsilon_hat^2
    
    #--------------------------------------------------------------------------
    # Create approximate posterior samples for transformation inference
    
    # Number of samples for approximating the posterior
    n_approx = 1000
    
    # Generate samples from approximate posterior at each observation
    # This mimics what pilot BART would provide
    z_samples = matrix(NA, n_approx, n)
    
    for(i in 1:n){
      # Sample from N(rf_mean[i], rf_var_total[i])
      z_samples[,i] = rnorm(n_approx, 
                            mean = rf_mean[i], 
                            sd = sqrt(rf_var_total[i]))
    }
    
    # Compute CDF approximations needed for transformation
    # Grid for z values
    z_grid = seq(min(z_samples) - 2*sigma_epsilon_hat, 
                 max(z_samples) + 2*sigma_epsilon_hat, 
                 length.out = ngrid)
    
    # Compute F_{Z|X=x_i} for each observation
    Fzx_eval = matrix(NA, ngrid, n)
    for(i in 1:n){
      Fzx_eval[,i] = ecdf(z_samples[,i])(z_grid)
    }
    
    # CDF of z using RF estimates (initial approximation):
    Fz_eval = rowMeans(Fzx_eval)
    
    if(verbose) cat('Random Forest approximation complete. Initializing main sampler...\n')
  }

  # Compute initial transformation:
  g = SeBR:::g_fun(y = y0,
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
  g_inv = SeBR:::g_inv_approx(g = g, t_grid = y_grid)

  #----------------------------------------------------------------------------
  # Set up main BART sampler with the latent response on the identified scale
  X_df <- as.data.frame(X)
  colnames(X_df) <- paste0("X", 1:ncol(X))
  main_data <- cbind(data.frame(z = z), X_df)

  main_sampler = dbarts::dbarts(
    formula = z ~ .,
    data = main_data,
    control = dbarts::dbartsControl(
      verbose = FALSE,
      keepTrainingFits = TRUE,
      keepTrees = var_select,  # Keep trees only if variable selection is enabled
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
  post_g_raw = array(NA, c(nsave, length(y0)))
  post_g_mu = numeric(nsave)
  post_g_sigma = numeric(nsave)
  
  # Initialize variable selection tracking arrays
  if(var_select) {
    var_usage_count = rep(0, p)  # Track total usage of each variable
    var_importance_scores = rep(0, p)  # Track importance scores
    selected_vars = 1:p  # Initially all variables selected
  }

  # Initialize timer for progress reporting:
  if(verbose) timer0 = proc.time()[3]

  #----------------------------------------------------------------------------
  # MCMC:
  for(nsi in 1:(nburn + nsave)){

    #--------------------------------------------------------------------------
    # Sample the transformation
    # Bayesian bootstrap for Y:
    Fy_eval = (n/(n + 1)) * SeBR:::bb(y)(y0)

    # CDF of Z:
    if(!fixedX){
      # Random X: re-randomize the weights
      weights_x = rgamma(n = n, shape = 1)
      weights_x = weights_x/sum(weights_x)
      Fz_eval = Fzx_eval %*% weights_x
    } else {
      # Fixed X: simple average
      Fz_eval = rowMeans(Fzx_eval)
    }

    # Sample the transformation:
    g = SeBR:::g_fun(y = y0,
              Fy_eval = Fy_eval,
              z = z_grid,
              Fz_eval = Fz_eval)

    # Apply the transformation:
    z = g(y)

    # Update the inverse transformation:
    g_inv = SeBR:::g_inv_approx(g = g, t_grid = y_grid)

    #--------------------------------------------------------------------------
    # Update BART given the current latent response
    main_sampler$setResponse(z)
    current_sample = main_sampler$run(numBurnIn = 0, numSamples = 1)

    #----------------------------------------------------------------------------
    # Variable selection: extract trees and count variable usage
    if(var_select && nsi > nburn/2) {  # Start variable selection after half burn-in
      tryCatch({
        # Extract trees from current BART fit
        trees = main_sampler$getTrees()
        
        # Count variable usage in all trees
        if(!is.null(trees) && length(trees) > 0) {
          for(tree_idx in 1:length(trees)) {
            tree = trees[[tree_idx]]
            if(!is.null(tree) && nrow(tree) > 0) {
              # Extract variable indices used in splits (excluding leaf nodes)
              split_vars = tree$var[tree$var > 0]  # var > 0 indicates split nodes
              if(length(split_vars) > 0) {
                # Count usage of each variable
                for(var_idx in split_vars) {
                  if(var_idx <= p) {  # Ensure valid variable index
                    var_usage_count[var_idx] = var_usage_count[var_idx] + 1
                  }
                }
              }
            }
          }
        }
      }, error = function(e) {
        # Silently handle any errors in tree extraction
        if(verbose && nsi %% 100 == 0) {
          cat("Warning: Could not extract trees for variable selection at iteration", nsi, "\n")
        }
      })
    }

    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Predictive samples:
      # Get BART predictions at test points using object method
      X_test_df <- as.data.frame(X_test)
      colnames(X_test_df) <- paste0("X", 1:ncol(X_test))
      f_hat = as.numeric(main_sampler$predict(X_test_df))

      sigma_bart = current_sample$sigma[1]
      ztilde = rnorm(n_test, mean = f_hat, sd = sigma_bart)
      
      post_ypred[nsi - nburn,] = g_inv(ztilde)

      g_values <- g(y0)
      train_fits <- current_sample$train
      mu_bart <- mean(as.numeric(train_fits))
      if (!is.finite(mu_bart)) {
        mu_bart <- 0
      }
      if (!is.finite(sigma_bart) || sigma_bart < 1e-8) {
        sigma_bart <- 1e-8
      }

      post_g_raw[nsi - nburn,] = g_values
      post_g_mu[nsi - nburn] = mu_bart
      post_g_sigma[nsi - nburn] = sigma_bart
      post_g[nsi - nburn,] = (g_values - mu_bart) / sigma_bart
    }

    if(verbose) SeBR:::computeTimeRemaining(nsi, timer0, nsave + nburn)
  }

  #----------------------------------------------------------------------------
  # Variable selection: compute importance scores and select variables
  if(var_select) {
    # Compute variable importance scores (normalized usage counts)
    total_usage = sum(var_usage_count)
    if(total_usage > 0) {
      var_importance_scores = var_usage_count / total_usage
    } else {
      var_importance_scores = rep(1/p, p)  # Equal importance if no usage tracked
    }
    
    # Select variables above threshold
    selected_vars = which(var_importance_scores > var_select_threshold)
    
    # Ensure at least one variable is selected
    if(length(selected_vars) == 0) {
      selected_vars = which.max(var_importance_scores)
    }
    
    if(verbose) {
      cat("\nVariable Selection Results:\n")
      cat("  Threshold:", var_select_threshold, "\n")
      cat("  Selected variables:", length(selected_vars), "of", p, "\n")
      cat("  Variable indices:", paste(selected_vars, collapse = ", "), "\n")
      cat("  Importance scores:", paste(round(var_importance_scores, 4), collapse = ", "), "\n")
    }
  }

  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer0
    if(tfinal > 60){
      cat(paste('Total time:', round(tfinal/60,1), 'minutes\n'))
    } else cat(paste('Total time:', round(tfinal), 'seconds\n'))
  }

  # Prepare return list
  result_list = list(
    fitted.values = colMeans(post_ypred),
    post_ypred = post_ypred,
    post_g = post_g,
    post_g_raw = post_g_raw,
    post_g_mu = post_g_mu,
    post_g_sigma = post_g_sigma,
    model = 'sebart', 
    y = y, X = X, X_test = X_test, 
    fixedX = fixedX,
    pilot_method = pilot_method,
    ntree = ntree,
    n.threads = n.threads, seed = seed,
    var_select = var_select
  )
  
  # Add variable selection results if enabled
  if(var_select) {
    result_list$var_importance_scores = var_importance_scores
    result_list$selected_vars = selected_vars
    result_list$var_select_threshold = var_select_threshold
  }
  
  # Add method-specific parameters to return list
  if(pilot_method == "bart") {
    result_list$pilot_ndraws = pilot_ndraws
    result_list$pilot_nburn = pilot_nburn
  } else {
    result_list$rf_fit = rf_fit
    result_list$rf_ntree = rf_ntree
    result_list$rf_nodesize = rf_nodesize
    result_list$rf_mtry = rf_mtry
  }
  
  return(result_list)
}
