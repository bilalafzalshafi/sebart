#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

read_dataset_raw <- function(dataset_name) {
  if (dataset_name == "forestfires") {
    read.csv("forestfires.csv")
  } else if (dataset_name == "bikesharing") {
    read.csv("hour.csv")
  } else if (dataset_name == "parkinsons") {
    data <- read.csv("parkinsons_updrs.data")
    data$motor_norm <- data$motor_UPDRS / 100
    data$total_norm <- data$total_UPDRS / 200
    data
  } else {
    stop("Unknown dataset: ", dataset_name)
  }
}

extract_target_values <- function(dataset_name, target_feature) {
  data <- read_dataset_raw(dataset_name)
  if (!target_feature %in% names(data)) {
    stop("Target feature ", target_feature, " not found in dataset ", dataset_name)
  }
  y_raw <- data[[target_feature]]
  y_raw <- y_raw[is.finite(y_raw) & !is.na(y_raw)]
  if (length(y_raw) == 0) {
    stop("No valid observations for ", dataset_name, " target ", target_feature)
  }
  data.frame(
    dataset = dataset_name,
    target = target_feature,
    y_raw = y_raw
  )
}

datasets_to_plot <- list(
  forestfires = c("area"),
  bikesharing = c("windspeed"),
  parkinsons = c("DFA")
)

all_targets <- do.call(rbind, lapply(names(datasets_to_plot), function(dataset_name) {
  targets <- datasets_to_plot[[dataset_name]]
  do.call(rbind, lapply(targets, function(target_feature) {
    extract_target_values(dataset_name, target_feature)
  }))
}))

all_targets$facet_label <- paste0(tools::toTitleCase(all_targets$dataset), " (", all_targets$target, ")")

summary_table <- all_targets %>%
  group_by(dataset, target) %>%
  summarise(
    n = n(),
    mean_raw = mean(y_raw),
    median_raw = median(y_raw),
    sd_raw = sd(y_raw),
    min_raw = min(y_raw),
    max_raw = max(y_raw),
    q1_raw = quantile(y_raw, 0.25),
    q3_raw = quantile(y_raw, 0.75),
    skew_raw = {
      sd_val <- sd(y_raw)
      if (sd_val > 0) mean(((y_raw - mean(y_raw)) / sd_val)^3) else NA_real_
    },
    .groups = "drop"
  )

write.csv(summary_table, file = "real_data_y_summary.csv", row.names = FALSE)

raw_plot <- ggplot(all_targets, aes(x = y_raw)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#f4b400", color = "white", alpha = 0.85) +
  geom_density(color = "#b45309", linewidth = 0.7) +
  facet_wrap(~ facet_label, ncol = 3, scales = "free") +
  labs(x = "Raw response", y = "Density",
       title = "Raw response marginals across real-data targets") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave("real_data_y_marginals.png", raw_plot,
       width = 9, height = 3.3, dpi = 300, bg = "white")

cat("Saved plot: real_data_y_marginals.png\n")
cat("Summary table: real_data_y_summary.csv\n")
