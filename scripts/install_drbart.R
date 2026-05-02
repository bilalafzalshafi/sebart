#!/usr/bin/env Rscript
# Installs DR-BART from a pinned commit so simulation results remain auditable.
# Usage: Rscript scripts/install_drbart.R

sha <- "e3c16129d52155b909b5189863aacb98c223799b"
repo <- "bilalafzalshafi/drbart"

cat(sprintf("Installing %s at commit %s...\n", repo, sha))
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(sprintf("%s@%s", repo, sha), upgrade = "never", dependencies = TRUE)
cat("Installation complete.\n")
