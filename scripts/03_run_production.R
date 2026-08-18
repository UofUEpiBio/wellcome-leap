source("config/simulation.R")
source("R/build_model.R")
source("R/extract_outcomes.R")
source("R/simulate.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1L) as.integer(args[[1]]) else 5000L
workers <- if (length(args) >= 2L) {
  as.integer(args[[2]])
} else {
  min(8L, parallel::detectCores())
}
batch_size <- if (length(args) >= 3L) as.integer(args[[3]]) else 100L

config <- default_simulation_config()
if (file.exists("data/derived/calibration.rds")) {
  config <- readRDS("data/derived/calibration.rds")$config
}

manifest <- run_simulation_batches(
  config,
  n_reps      = n_reps,
  batch_size  = batch_size,
  workers     = workers,
  output_dir  = "data/derived/production"
)
print(manifest)
