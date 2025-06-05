# basic_bart.R - experimenting with BART

library(dbarts)
library(BayesTree)
library(BART)

set.seed(42)

# different function types
n <- 200
x <- seq(-2, 2, length.out = n)
X <- matrix(x, ncol = 1)

# linear function
y_linear <- 2 * x + 1 + rnorm(n, 0, 0.3)
bart_linear <- bart(X, y_linear, ndpost = 500, nskip = 100, verbose = FALSE)
pred_linear <- colMeans(bart_linear$yhat.train)
rmse_linear <- sqrt(mean((pred_linear - (2 * x + 1))^2))
cat("Linear RMSE:", rmse_linear, "\n")

# sine function
y_sine <- 3 * sin(2 * pi * x) + rnorm(n, 0, 0.3)
bart_sine <- bart(X, y_sine, ndpost = 500, nskip = 100, verbose = FALSE)
pred_sine <- colMeans(bart_sine$yhat.train)
rmse_sine <- sqrt(mean((pred_sine - 3 * sin(2 * pi * x))^2))
cat("Sine RMSE:", rmse_sine, "\n")

# step function
y_step <- ifelse(x > 0, 2, -1) + rnorm(n, 0, 0.3)
bart_step <- bart(X, y_step, ndpost = 500, nskip = 100, verbose = FALSE)
pred_step <- colMeans(bart_step$yhat.train)
rmse_step <- sqrt(mean((pred_step - ifelse(x > 0, 2, -1))^2))
cat("Step RMSE:", rmse_step, "\n")

# check uncertainty
y_test <- sin(2 * x) + rnorm(n, 0, 0.4)
bart_uncertain <- bart(X, y_test, ndpost = 800, nskip = 200, verbose = FALSE)
pred_mean <- colMeans(bart_uncertain$yhat.train)
pred_lower <- apply(bart_uncertain$yhat.train, 2, quantile, 0.025)
pred_upper <- apply(bart_uncertain$yhat.train, 2, quantile, 0.975)
coverage <- mean((y_test >= pred_lower) & (y_test <= pred_upper))
cat("Coverage:", coverage, "\n")

# train/test split
n_split <- 300
x_all <- runif(n_split, -1, 1)
y_true <- x_all^2 + sin(pi * x_all)
y_all <- y_true + rnorm(n_split, 0, 0.3)

train_idx <- sample(n_split, 200)
test_idx <- setdiff(1:n_split, train_idx)

X_train <- matrix(x_all[train_idx], ncol = 1)
X_test <- matrix(x_all[test_idx], ncol = 1)
y_train <- y_all[train_idx]
y_test <- y_all[test_idx]

bart_split <- bart(X_train, y_train, x.test = X_test, ndpost = 600, nskip = 150, verbose = FALSE)
pred_test <- colMeans(bart_split$yhat.test)
test_rmse <- sqrt(mean((pred_test - (x_all[test_idx]^2 + sin(pi * x_all[test_idx])))^2))
cat("Test RMSE:", test_rmse, "\n")

# multidimensional
n_multi <- 250
x1 <- runif(n_multi, -1, 1)
x2 <- runif(n_multi, -1, 1)
X_multi <- cbind(x1, x2)
y_multi <- x1^2 + x2 + x1*x2 + rnorm(n_multi, 0, 0.2)

bart_multi <- bart(X_multi, y_multi, ndpost = 500, nskip = 100, verbose = FALSE)
pred_multi <- colMeans(bart_multi$yhat.train)
multi_rmse <- sqrt(mean((pred_multi - (x1^2 + x2 + x1*x2))^2))
cat("Multi RMSE:", multi_rmse, "\n")

save(bart_linear, bart_sine, bart_step, bart_multi, file = "bart_experiments.RData")