source("config/multiscale.R")
source("R/within_host.R")
source("R/omics.R")
source("R/reproduction.R")
source("R/multiscale.R")

args <- commandArgs(trailingOnly = TRUE)
profiles_per_site <- if (length(args)) as.integer(args[[1]]) else 4L
cores <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L
config <- default_multiscale_config()
study <- simulate_multiscale_dataset(
  config,
  profiles_per_site = profiles_per_site,
  fit_profiles = FALSE,
  cores = cores
)
dir.create("data/derived/multiscale", recursive = TRUE, showWarnings = FALSE)
saveRDS(study, "data/derived/multiscale/truth_study.rds", compress = "gzip")
targets <- combine_multiscale_targets(study$truth)
print(utils::head(targets, 12L))
cat(
  "Saved",
  length(study$simulated),
  "profiles and",
  nrow(targets),
  "profile/scenario target rows.\n"
)
