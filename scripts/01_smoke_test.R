source("config/simulation.R")
source("R/build_model.R")
source("R/extract_outcomes.R")
source("R/simulate.R")

config <- default_simulation_config()
result <- run_simulation_study(config, n_reps = 2L, workers = 4L)

print(result$runs)
print(aggregate(
  secondary_cases ~ scenario,
  filter_analysis_agents(result$agents),
  mean
))

