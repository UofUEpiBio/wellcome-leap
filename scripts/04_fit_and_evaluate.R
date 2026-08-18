source("config/simulation.R")
source("R/build_model.R")
source("R/extract_outcomes.R")
source("R/simulate.R")
source("R/masking.R")
source("R/torch_model.R")
source("R/fit_models.R")

manifest_path <- "data/derived/manifest.rds"
if (!file.exists(manifest_path)) {
  stop("Run scripts/03_run_production.R before fitting the ML model.")
}
manifest <- readRDS(manifest_path)
batches <- lapply(manifest$path, readRDS)
agents <- do.call(rbind, lapply(batches, `[[`, "agents"))

fit <- fit_masked_model(agents)
mae <- evaluate_masked_mae(fit)
save_masked_model(fit)
saveRDS(mae, "artifacts/masked_model_mae.rds")
print(fit$history)
print(mae)
