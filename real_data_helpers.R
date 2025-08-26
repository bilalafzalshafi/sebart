library(datasets)
if (!require(moments, quietly = TRUE)) {
  install.packages("moments")
  library(moments)
}

# Load and preprocess real datasets
# @param dataset_name Name of dataset 
# @param target_feature Which feature to use as response variable
# @param n_train Number of training observations
# @param n_test Number of test observations  
# @param seed Random seed for reproducibility
load_real_dataset <- function(dataset_name, target_feature = NULL, 
                             n_train = 200, n_test = 500, seed = 123) {
  
  set.seed(seed)
  
  if (dataset_name == "forestfires") {
    data <- read.csv("forestfires.csv")
    
    if (is.null(target_feature)) target_feature <- "area"
    y_raw <- data[[target_feature]]
    X_raw <- data[, !names(data) %in% c(target_feature)]

  } 
  else if (dataset_name == "bikesharing") {
    data <- read.csv("hour.csv")

    if (is.null(target_feature)) target_feature <- "hum"  # Changed default from "temp" to "hum"
    y_raw <- data[[target_feature]]
    
    exclude_cols <- c(target_feature, "instant", "dteday", "yr", "casual", "registered", "cnt")
    if (target_feature != "temp") exclude_cols <- c(exclude_cols, "temp")
    if (target_feature != "hum") exclude_cols <- c(exclude_cols, "hum") 
    if (target_feature != "windspeed") exclude_cols <- c(exclude_cols, "windspeed")
    
    X_raw <- data[, !names(data) %in% exclude_cols]
  }

  else if (dataset_name == "parkinsons") {
    data <- read.csv("parkinsons_updrs.data")
    
    data$motor_norm <- data$motor_UPDRS / 100
    data$total_norm <- data$total_UPDRS / 200
    
    if (is.null(target_feature)) target_feature <- "NHR"
    y_raw <- data[[target_feature]]
    
    exclude_cols <- c(target_feature, "motor_UPDRS", "total_UPDRS", 
                      "motor_norm", "total_norm")
    X_raw <- data[, !names(data) %in% exclude_cols]

  }

  else {
    stop("Unknown dataset: ", dataset_name, 
         ". Available: ", paste(get_available_real_datasets(), collapse = ", "))
  }
  
  X_raw <- X_raw[, sapply(X_raw, is.numeric), drop = FALSE]
  
  if (any(is.na(y_raw))) {
    complete_idx <- !is.na(y_raw) & complete.cases(X_raw)
    y_raw <- y_raw[complete_idx]
    X_raw <- X_raw[complete_idx, , drop = FALSE]
  }
  
  # Normalize y values to [0, 1] range
  y_min <- min(y_raw)
  y_max <- max(y_raw)
  y_normalized <- (y_raw - y_min) / (y_max - y_min)
  
  X_scaled <- scale(X_raw)
  
  # Train/test split
  n_total <- length(y_normalized)
  available_train <- min(n_train, floor(n_total * 0.7))
  available_test <- min(n_test, n_total - available_train)
  
  train_idx <- sample(1:n_total, available_train)
  remaining_idx <- setdiff(1:n_total, train_idx)
  test_idx <- sample(remaining_idx, available_test)
  
  return(list(
    y_train = y_normalized[train_idx],
    y_test = y_normalized[test_idx],
    X_train = X_scaled[train_idx, , drop = FALSE],
    X_test = X_scaled[test_idx, , drop = FALSE],
    f_true_train = NULL,
    f_true_test = NULL,
    z_train = y_normalized[train_idx],
    z_test = y_normalized[test_idx],
    g_true = NULL,
    y_unique = sort(unique(y_normalized)),
    dataset_info = list(
      name = dataset_name,
      target = target_feature,
      n_features = ncol(X_scaled),
      n_total = n_total,
      y_range = c(y_min, y_max)
    )
  ))
}

get_available_real_datasets <- function() {
  c("forestfires", "bikesharing", "parkinsons")
}

get_target_features <- function(dataset_name) {
  switch(dataset_name,
    "forestfires" = c("area"),
    "bikesharing" = c("hum", "windspeed"),  # Removed "temp"
    "parkinsons" = c("NHR", "HNR", "DFA", "motor_norm", "total_norm"),
    stop("Unknown dataset: ", dataset_name)
  )
}

demo_real_data_diagnostics <- function(dataset_name, target_feature = NULL, 
                                      n_train = 200, n_test = 500) {
  
  real_data <- load_real_dataset(dataset_name, target_feature, 
                                n_train, n_test, seed = 123)
  
  # Updated to use sebart instead of sebart
  sebart_fit <- sebart(y = real_data$y_train, 
                     X = real_data$X_train, 
                     X_test = real_data$X_test,
                     ntree = 200, nsave = 1000, nburn = 1000, 
                     verbose = FALSE)
  
  metrics <- create_all_sebart_diagnostics(sebart_fit, 
                                         real_data$y_test, 
                                         true_g_values = NULL)
  
  return(list(
    sebart_fit = sebart_fit,
    real_data = real_data,
    metrics = metrics,
    dataset_info = real_data$dataset_info
  ))
}

summarize_target_distribution <- function(dataset_name, target_feature = NULL) {
  
  real_data <- load_real_dataset(dataset_name, target_feature, 
                                n_train = 100, n_test = 100, seed = 42)
  
  y_combined <- c(real_data$y_train, real_data$y_test)
  
  summary_stats <- list(
    dataset = dataset_name,
    target = target_feature,
    n_obs = length(y_combined),
    mean = mean(y_combined),
    median = median(y_combined),
    sd = sd(y_combined),
    skewness = moments::skewness(y_combined),
    min = min(y_combined),
    max = max(y_combined),
    q25 = quantile(y_combined, 0.25),
    q75 = quantile(y_combined, 0.75)
  )
  
  return(summary_stats)
}