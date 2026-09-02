source("config/multiscale.R")
source("R/within_host.R")
source("R/omics.R")
source("R/reproduction.R")
source("R/multiscale.R")

args <- commandArgs(trailingOnly = TRUE)
profiles_per_site <- if (length(args)) as.integer(args[[1]]) else 4L
config <- default_multiscale_config()
study <- simulate_multiscale_dataset(
  config,
  profiles_per_site = profiles_per_site,
  fit_profiles = FALSE
)
dir.create("data/derived/multiscale", recursive = TRUE, showWarnings = FALSE)
saveRDS(study, "data/derived/multiscale/truth_study.rds", compress = "xz")
print(combine_multiscale_targets(study$truth))
