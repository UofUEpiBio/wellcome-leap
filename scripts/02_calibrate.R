source("config/simulation.R")
source("R/build_model.R")
source("R/extract_outcomes.R")
source("R/simulate.R")
source("R/calibrate.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1L) as.integer(args[[1]]) else 200L
workers <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L

config <- default_simulation_config()
calibration <- calibrate_scenarios(config, n_reps = n_reps, workers = workers)
dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
saveRDS(calibration, "data/derived/calibration.rds")
print(calibration$history)
print(calibration$config[c("alpha", "delta_scenario")])

