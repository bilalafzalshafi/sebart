#!/usr/bin/env Rscript

resolve_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) && nzchar(file_arg)) {
    script_dir <- dirname(normalizePath(file_arg))
  } else {
    ofile <- NULL
    frame1 <- sys.frames()[[1]]
    if (!is.null(frame1$ofile)) {
      ofile <- frame1$ofile
    } else if (!is.null(frame1$`__file__`)) {
      ofile <- frame1$`__file__`
    }
    script_dir <- if (!is.null(ofile)) dirname(normalizePath(ofile)) else normalizePath(getwd())
  }
  if (basename(script_dir) %in% c("scripts", "bin")) {
    normalizePath(file.path(script_dir, ".."))
  } else {
    script_dir
  }
}

parse_args <- function(args) {
  parsed <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    value <- if (length(kv) > 1) kv[2] else NA_character_
    parsed[[key]] <- value
  }
  parsed
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)) a else b

as_bool <- function(x) tolower(x) %in% c("true", "t", "1", "yes", "y")

format_seconds <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 60) return(sprintf("%.2fs", x))
  if (x < 3600) return(sprintf("%.2fmin", x/60))
  sprintf("%.2fh", x/3600)
}

main <- function() {
  root <- resolve_root()
  options(sebart.root = root)
  source(file.path(root, "R", "experiment_utils.R"))
  source(file.path(root, "sebart.R"))
  source(file.path(root, "bart_bc.R"))

  defaults <- list(
    scenario = "sigmoid",
    n_train = "200",
    n_test = "1000",
    p = "10",
    models = "sblm,bart_bc,bart,sebart",
    seed = "202401",
    predictor = "gaussian",
    alpha = "0.1"
  )

  args <- modifyList(defaults, parse_args(commandArgs(trailingOnly = TRUE)))
  scenario <- args$scenario
  n_train <- as.integer(args$n_train)
  n_test <- as.integer(args$n_test)
  p <- as.integer(args$p)
  models <- trimws(strsplit(args$models, ",")[[1]])
  seed <- as.integer(args$seed)
  predictor <- args$predictor
  alpha <- as.numeric(args$alpha)

  config <- list(
    seed = seed,
    predictor = predictor,
    model = list(ntree = 200, nsave = 1000, nburn = 1000, verbose = FALSE)
  )

  iter_res <- run_simulation_iteration(
    scenario = scenario,
    iteration = 1,
    total_iterations = 1,
    config = config,
    models = models,
    n_train = n_train,
    n_test = n_test,
    p = p,
    seed_offset = 0,
    alpha = alpha,
    record_transformation = FALSE,
    verbose_progress = TRUE
  )

  metrics <- iter_res$metrics
  metrics <- metrics[metrics$Model %in% models, , drop = FALSE]
  if (!nrow(metrics)) {
    stop("No metrics returned (check models list).", call. = FALSE)
  }
  metrics$RuntimeFmt <- vapply(metrics$Runtime, format_seconds, character(1))
  metrics <- metrics[match(models, metrics$Model), , drop = FALSE]

  cat(sprintf("Scenario: %s, n_train=%d, n_test=%d, seed=%d\n", scenario, n_train, n_test, seed))
  cat("Per-model training + prediction times:\n")
  for (idx in seq_len(nrow(metrics))) {
    cat(sprintf("  %-8s : %s\n", metrics$Model[idx], metrics$RuntimeFmt[idx]))
  }
}

if (sys.nframe() == 0) {
  main()
}
