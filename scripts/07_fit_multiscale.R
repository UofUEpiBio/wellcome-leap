source("config/multiscale.R")
source("R/within_host.R")
source("R/omics.R")
source("R/reproduction.R")
source("R/multiscale.R")
source("R/multiscale_abm.R")

truth_path <- "data/derived/multiscale/truth_study.rds"
if (!file.exists(truth_path)) {
  stop("Run scripts/06_simulate_multiscale.R before fitting profiles.")
}
args <- commandArgs(trailingOnly = TRUE)
n_starts <- if (length(args)) as.integer(args[[1]]) else 3L
config <- default_multiscale_config()
study <- readRDS(truth_path)
study$fitted <- lapply(seq_along(study$simulated), function(index) {
  simulated <- study$simulated[[index]]
  fit <- fit_within_host_profile(
    simulated$observations,
    config,
    antibiotic_end = simulated$profile$antibiotic_end,
    n_starts = n_starts,
    seed = config$base_seed + 10000L + index
  )
  build_multiscale_profile_bundle(
    simulated,
    config,
    parameter_source = "fitted",
    fit = fit
  )
})

profile_rows <- lapply(study$fitted, function(bundle) {
  summary <- bundle$baseline$between$re
  data.frame(
    profile_id = bundle$profile$profile_id,
    effective_contact_risk = summary[["effective_contact_risk"]],
    carriage_duration = summary[["carriage_duration"]],
    contact_rate = bundle$profile$contact_rate,
    susceptible_fraction = bundle$profile$susceptible_fraction
  )
})
profile_data <- do.call(rbind, profile_rows)
abm <- run_multiscale_abm(
  profile_data,
  n_agents = 2000L,
  prevalence = 5 / 2000,
  days = 120L,
  seed = config$base_seed + 20000L
)
study$abm <- abm
saveRDS(study, "data/derived/multiscale/fitted_study.rds", compress = "xz")
print(combine_multiscale_targets(study$fitted))
