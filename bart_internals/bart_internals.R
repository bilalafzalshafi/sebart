library(dbarts)
library(BART)

set.seed(42)

# basic test data
n <- 80
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
X <- cbind(x1, x2)
y <- 2 * sin(pi * x1) + x2^2 + rnorm(n, 0, 0.3)

# compare regular BART vs dbarts creation
bart_regular <- bart(X, y, ndpost = 200, nskip = 50, verbose = FALSE)
cat("Regular BART components:", names(bart_regular), "\n")
cat("Regular BART predictions dim:", dim(bart_regular$yhat.train), "\n")

# create dbarts sampler
sampler <- dbarts(y ~ x1 + x2, data = data.frame(y = y, x1 = x1, x2 = x2))
cat("dbarts sampler components:", names(sampler), "\n")

# examine control structure
control_obj <- sampler$control
cat("Available control slots:", slotNames(control_obj), "\n")
cat("Number of trees:", control_obj@n.trees, "\n")
cat("Number of samples:", control_obj@n.samples, "\n")
cat("Burn-in:", control_obj@n.burn, "\n")

# examine data structure
data_obj <- sampler$data
cat("Data object class:", class(data_obj), "\n")
cat("Response variable length:", length(data_obj@y), "\n")
cat("Predictor matrix dim:", dim(data_obj@x), "\n")

# test different control configurations
config_small <- dbartsControl(n.trees = 25L, n.samples = 100L, n.burn = 50L)
config_large <- dbartsControl(n.trees = 300L, n.samples = 100L, n.burn = 50L)

sampler_small <- dbarts(y ~ x1 + x2, data = data.frame(y = y, x1 = x1, x2 = x2), 
                       control = config_small)
sampler_large <- dbarts(y ~ x1 + x2, data = data.frame(y = y, x1 = x1, x2 = x2), 
                       control = config_large)

cat("Small config trees:", sampler_small$control@n.trees, "\n")
cat("Large config trees:", sampler_large$control@n.trees, "\n")

# test basic run functionality
result_regular <- sampler$run()
cat("dbarts run result type:", class(result_regular), "\n")
if (is.list(result_regular)) {
  cat("dbarts result components:", names(result_regular), "\n")
}

# compare prediction approaches
pred_bart_regular <- colMeans(bart_regular$yhat.train)
rmse_bart_regular <- sqrt(mean((pred_bart_regular - y)^2))
cat("Regular BART RMSE:", rmse_bart_regular, "\n")

# test memory footprint
size_bart <- object.size(bart_regular)
size_sampler <- object.size(sampler)
cat("Regular BART size:", format(size_bart, units = "Kb"), "\n")
cat("dbarts sampler size:", format(size_sampler, units = "Kb"), "\n")

# test parameter access consistency
configs_test <- list(
  tiny = dbartsControl(n.trees = 10L),
  medium = dbartsControl(n.trees = 100L),
  large = dbartsControl(n.trees = 500L)
)

for (config_name in names(configs_test)) {
  test_sampler <- dbarts(y ~ x1 + x2, data = data.frame(y = y, x1 = x1, x2 = x2),
                        control = configs_test[[config_name]])
  n_trees <- test_sampler$control@n.trees
  cat(config_name, "configuration trees:", n_trees, "\n")
}

# test data modification workflow
y_modified <- y + rnorm(n, 0, 0.1)
data_modified <- data.frame(y_modified = y_modified, x1 = x1, x2 = x2)

sampler_modified <- dbarts(y_modified ~ x1 + x2, data = data_modified)
cat("Modified data sampler created successfully\n")
cat("Modified data response length:", length(sampler_modified$data@y), "\n")

# compare with transformed data
y_transformed <- log(abs(y) + 1)  # simple transformation
data_transformed <- data.frame(y_transformed = y_transformed, x1 = x1, x2 = x2)

sampler_transformed <- dbarts(y_transformed ~ x1 + x2, data = data_transformed)
cat("Transformed data sampler created successfully\n")

# basic capability assessment
capabilities <- list(
  sampler_creation = TRUE,
  control_parameter_access = !is.null(sampler$control@n.trees),
  data_object_access = !is.null(sampler$data@y),
  multiple_configurations = length(configs_test) > 0,
  data_modification = !is.null(sampler_modified),
  transformation_workflow = !is.null(sampler_transformed)
)

cat("dbarts capabilities:\n")
for (cap_name in names(capabilities)) {
  cat("  ", cap_name, ":", capabilities[[cap_name]], "\n")
}

save(sampler, bart_regular, capabilities, configs_test, 
     file = "bart_internals.RData")