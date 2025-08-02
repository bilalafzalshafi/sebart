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

    if (is.null(target_feature)) target_feature <- "temp"
    y_raw <- data[[target_feature]]
    
    exclude_cols <- c(target_feature, "instant", "dteday", "yr", "casual", "registered", "cnt")
    if (target_feature != "temp") exclude_cols <- c(exclude_cols, "temp")
    if (target_feature != "hum") exclude_cols <- c(exclude_cols, "hum") 
    if (target_feature != "windspeed") exclude_cols <- c(exclude_cols, "windspeed")
    
    X_raw <- data[, !names(data) %in% exclude_cols]
  }

  else if (dataset_name == "concrete") {
    data <- read.csv("Concrete_Data.csv")

    names(data) <- c("cement", "slag", "flyash", "water", "superplast", 
                "coarse_agg", "fine_agg", "age", "strength")
    
    if (is.null(target_feature)) target_feature <- "strength"
    y_raw <- data[[target_feature]]
    X_raw <- data[, !names(data) %in% c(target_feature)]
  }

  else if (dataset_name == "studentperformance") {
    data <- read.csv("student-mat.csv", sep = ";")
    
    data$G3_norm <- data$G3 / 20  # Normalize to [0,1]
    data$absences_norm <- data$absences / max(data$absences)
    data$failures_norm <- data$failures / max(data$failures)
    
    if (is.null(target_feature)) target_feature <- "G3_norm"
    y_raw <- data[[target_feature]]
    
    exclude_cols <- c(target_feature, "G1", "G2", "G3", "G3_norm", 
                      "absences_norm", "failures_norm")
    X_raw <- data[, !names(data) %in% exclude_cols]

  }

  else if (dataset_name == "glass") {
      data <- read.csv("glass.data", header = FALSE)
      names(data) <- c("Id", "RI", "Na", "Mg", "Al", "Si", "K", "Ca", "Ba", "Fe", "Type")
      
      if (is.null(target_feature)) target_feature <- "Na"
      
      oxide_cols <- c("Na", "Mg", "Al", "Si", "K", "Ca", "Ba", "Fe")
      for (col in oxide_cols) {
        data[[col]] <- data[[col]] / 100
      }
      
      y_raw <- data[[target_feature]]
      exclude_cols <- c("Id", target_feature, "Type")
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

  else if (dataset_name == "yeast") {
    data <- read.table("yeast.data")
    names(data) <- c("seq_name", "mcg", "gvh", "alm", "mit", "erl", 
                     "pox", "vac", "nuc", "class")
    
    if (is.null(target_feature)) target_feature <- "mcg"
    y_raw <- data[[target_feature]]
    
    exclude_cols <- c("seq_name", target_feature, "class")
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
      n_total = n_total
    )
  ))
}

get_available_real_datasets <- function() {
  c("forestfires", "bikesharing", "concrete", "studentperformance", "glass", 
  "parkinsons", "yeast")
}

get_target_features <- function(dataset_name) {
  switch(dataset_name,
    "forestfires" = c("area"),
    "bikesharing" = c("temp", "hum", "windspeed"),
    "concrete" = c("strength"),
    "studentperformance" = c("G3_norm", "absences_norm", "failures_norm"),
    "glass" = c("Na", "Mg", "Al", "Si", "K", "Ca", "Ba", "Fe"),
    "parkinsons" = c("NHR", "HNR", "DFA", "motor_norm", "total_norm"),
    "yeast" = c("mcg", "gvh", "alm", "mit", "pox", "vac", "nuc"),
    stop("Unknown dataset: ", dataset_name)
  )
}

demo_real_data_diagnostics <- function(dataset_name, target_feature = NULL, 
                                      n_train = 200, n_test = 500) {
  
  real_data <- load_real_dataset(dataset_name, target_feature, 
                                n_train, n_test, seed = 123)
  
  sbart_fit <- sbart(y = real_data$y_train, 
                     X = real_data$X_train, 
                     X_test = real_data$X_test,
                     ntree = 200, nsave = 1000, nburn = 1000, 
                     verbose = FALSE)
  
  metrics <- create_all_sbart_diagnostics(sbart_fit, 
                                         real_data$y_test, 
                                         true_g_values = NULL)
  
  return(list(
    sbart_fit = sbart_fit,
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