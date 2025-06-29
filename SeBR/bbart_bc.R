#' Bayesian BART with a Box-Cox transformation
#'
#' MCMC sampling for Bayesian BART regression with a
#' (known or unknown) Box-Cox transformation.
#'
#' @param y \code{n x 1} response vector
#' @param X \code{n x p} matrix of predictors
#' @param X_test \code{n_test x p} matrix of predictors for test data;
#' default is the observed covariates \code{X}
#' @param ntree number of trees in BART ensemble
#' @param lambda Box-Cox transformation; if NULL, estimate this parameter
#' @param sample_lambda logical; if TRUE, sample lambda, otherwise
#' use the fixed value of lambda above or the MLE (if lambda unspecified)
#' @param nsave number of MCMC iterations to save
#' @param nburn number of MCMC iterations to discard
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' \item \code{post_lambda} \code{nsave} posterior samples of lambda
#' \item \code{model}: the model fit (here, \code{bbart_bc})
#' }
#' as well as the arguments passed in.
#'
#' @details This function provides Bayesian inference for
#' transformed BART models via MCMC sampling. The transformation is
#' parametric from the Box-Cox family, which has one parameter \code{lambda}.
#' That parameter may be fixed in advance or learned from the data.
#' For computational efficiency, the BART hyperparameters use standard
#' defaults and the algorithm alternates between updating the Box-Cox
#' parameter and the BART model.
#'
#' @note Box-Cox transformations may be useful in some cases, but
#' in general we recommend the nonparametric transformation (with
#' Monte Carlo sampling) in \code{\link{sbart}}.
#'
#' @examples
#' \donttest{
#' # Simulate some data:
#' n <- 200
#' p <- 5
#' X <- matrix(rnorm(n*p), n, p)
#' f_true <- sin(2*X[,1]) + X[,2]^2 - X[,3]*X[,4]
#' y <- g_inv_bc(f_true + rnorm(n, sd=0.1), lambda = 0.5) # square-root transformation
#' 
#' # Test data
#' X_test <- matrix(rnorm(50*p), 50, p)
#' 
#' # Fit Bayesian BART with Box-Cox transformation:
#' fit <- bbart_bc(y = y, X = X, X_test = X_test)
#' names(fit) # what is returned
#' round(quantile(fit$post_lambda), 3) # summary of unknown Box-Cox parameter
#' 
#' # Plot the model predictions (point and interval estimates):
#' pi_y <- t(apply(fit$post_ypred, 2, quantile, c(0.05, .95))) # 90% PI
#' plot(1:nrow(pi_y), pi_y[,1], type='n', ylim=range(pi_y),
#'      xlab='Test observation', ylab='y', main='Prediction intervals')
#' segments(1:nrow(pi_y), pi_y[,1], 1:nrow(pi_y), pi_y[,2])
#' points(1:nrow(pi_y), fitted(fit), pch=16)
#' }
#'
#' @export
bbart_bc = function(y, X, X_test = X,
                    ntree = 200,
                    lambda = NULL,
                    sample_lambda = TRUE,
                    nsave = 1000,
                    nburn = 1000,
                    nskip = 0,
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
  if(ntree < 10) stop('ntree must be at least 10')

  # Recurring terms:
  y0 = sort(unique(y))

  #----------------------------------------------------------------------------
  # Initialize the parameters:

  # Box-Cox parameter, if unspecified
  if(is.null(lambda)){

    # Conditional MLE using initial BART fit:
    # First get initial BART fit
    # Properly construct data frame with matrix columns
    X_df <- as.data.frame(X)
    colnames(X_df) <- paste0("X", 1:ncol(X))
    initial_data <- cbind(data.frame(y = y), X_df)
    
    # Create initial sampler (use defaults for priors)
    initial_sampler = dbarts::dbarts(
      formula = y ~ .,
      data = initial_data,
      control = dbarts::dbartsControl(
        verbose = FALSE,
        keepTrainingFits = TRUE,
        n.trees = ntree,
        n.chains = 1L,
        n.threads = 1L,
        rngSeed = seed
      )
    )

    initial_fit = initial_sampler$run(numBurnIn = 50, numSamples = 50)
    y_hat_initial = rowMeans(initial_fit$train)

    # MLE for lambda given initial fit:
    opt = optim(par = 0.5,
                fn = function(l_bc){
                  z_bc = g_bc(y, lambda = l_bc)
                  mean((z_bc - y_hat_initial)^2)
                }, method = "L-BFGS-B", lower = 0.01, upper = 1.99
    )
    if(opt$convergence == 0){
      lambda = opt$par
    } else {
      warning('Optimization failed...setting lambda = 1/2')
      lambda = 1/2
    }
  }

  # Latent data:
  z = g_bc(y, lambda = lambda)
  z_scaled = (z - mean(z))/sd(z) # standardize for BART
  z_scale_mean = mean(z)
  z_scale_sd = sd(z)

  #----------------------------------------------------------------------------
  # Set up BART sampler
  # Properly construct data frame with matrix columns
  X_df <- as.data.frame(X)
  colnames(X_df) <- paste0("X", 1:ncol(X))
  bart_data <- cbind(data.frame(z_scaled = z_scaled), X_df)
  
  # Create main BART sampler (use defaults for priors)
  bart_sampler = dbarts::dbarts(
    formula = z_scaled ~ .,
    data = bart_data,
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

  # Initialize BART sampler using object method
  invisible(bart_sampler$run(numBurnIn = 1, numSamples = 1))

  #----------------------------------------------------------------------------
  # Store MCMC output:
  post_ypred = array(NA, c(nsave, n_test))
  post_g = array(NA, c(nsave, length(y0)))
  post_lambda = rep(NA, nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    #----------------------------------------------------------------------------
    # Block 1: sample the transformation
    if(sample_lambda){
      
      # Get current BART fit for likelihood evaluation
      current_sample = bart_sampler$run(numBurnIn = 0, numSamples = 1)
      current_fit = as.vector(current_sample$train) * z_scale_sd + z_scale_mean
      current_sigma = current_sample$sigma[1] * z_scale_sd

      # Sample lambda using slice sampling:
      lambda = uni.slice(x0 = lambda,
                         g = function(l_bc){
                           # Likelihood for lambda:
                           z_bc = g_bc(y, lambda = l_bc)
                           sum(dnorm(z_bc,
                                     mean = current_fit,
                                     sd = current_sigma, log = TRUE)) +
                             # Prior on lambda, truncated to [0, 2]
                             dnorm(l_bc, mean = 1/2, sd = 1/2, log = TRUE)
                         },
                         w = 1/2, m = 50, lower = 0, upper = 2)

      # Update z and rescale:
      z = g_bc(y, lambda = lambda)
      z_scale_mean = mean(z)
      z_scale_sd = sd(z)
      z_scaled = (z - z_scale_mean)/z_scale_sd

      # Update BART sampler with new response using object method:
      bart_sampler$setResponse(z_scaled)
    }

    #----------------------------------------------------------------------------
    # Block 2: sample BART model
    
    # Single MCMC step for BART using object method
    bart_sample = bart_sampler$run(numBurnIn = 0, numSamples = 1)

    #----------------------------------------------------------------------------
    # Store the MCMC:
    if(nsi > nburn){

      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Predictive samples of ytilde:
        # Get BART predictions at test points using object method
        X_test_df <- as.data.frame(X_test)
        colnames(X_test_df) <- paste0("X", 1:ncol(X_test))
        bart_pred_scaled = bart_sampler$predict(X_test_df)
        bart_pred = as.vector(bart_pred_scaled) * z_scale_sd + z_scale_mean
        
        # Add noise
        current_sigma_scaled = bart_sample$sigma[1]
        ztilde = bart_pred + current_sigma_scaled * z_scale_sd * rnorm(n = n_test)
        post_ypred[isave,] = g_inv_bc(ztilde, lambda = lambda)

        # Posterior samples of the transformation:
        post_g[isave,] = g_bc(y0, lambda = lambda)
        post_lambda[isave] = lambda

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot)
  }
  # Summarize computing time:
  if(verbose){
    tfinal = proc.time()[3] - timer0
    if(tfinal > 60){
      print(paste('Total time:', round(tfinal/60,1), 'minutes'))
    } else print(paste('Total time:', round(tfinal), 'seconds'))
  }

  return(list(
    fitted.values = colMeans(post_ypred),
    post_ypred = post_ypred,
    post_g = post_g, 
    post_lambda = post_lambda,
    model = 'bbart_bc', 
    y = y, X = X, X_test = X_test, 
    ntree = ntree,
    sample_lambda = sample_lambda,
    n.threads = n.threads, seed = seed))
}