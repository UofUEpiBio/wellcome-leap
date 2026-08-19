source("config/simulation.R")
source("R/build_model.R")
source("R/extract_outcomes.R")
source("R/simulate.R")
source("R/masking.R")
source("R/torch_model.R")
source("R/fit_models.R")

manifest_path <- "data/derived/production/manifest.rds"
if (!file.exists(manifest_path)) {
  stop("Run scripts/03_run_production.R before fitting the ML model.")
}
torch::torch_set_num_threads(min(8L, parallel::detectCores()))
study <- load_simulation_batches(manifest_path)
agents <- filter_analysis_agents(study$agents)
rm(study)
gc()

fit <- fit_masked_model(
  agents,
  hidden_dim_1 = 32L,
  hidden_dim_2 = 16L,
  batch_size   = 2048L,
  max_epochs   = 60L,
  patience     = 15L
)
rmse <- evaluate_masked_rmse(fit)
save_masked_model(fit)
saveRDS(rmse, "artifacts/masked_model_rmse.rds")
print(fit$history)
print(rmse)
