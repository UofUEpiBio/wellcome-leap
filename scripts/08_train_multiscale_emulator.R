source("R/multimodal_emulator.R")
source("R/multiscale.R")

study_path <- "data/derived/multiscale/fitted_study.rds"
if (!file.exists(study_path)) {
  stop("Run scripts/07_fit_multiscale.R before training the emulator.")
}
study <- readRDS(study_path)
data <- make_multiscale_emulator_data(study$fitted)
fit <- fit_multiscale_emulator(data)
evaluation <- evaluate_multiscale_emulator(fit)
save_multiscale_emulator(fit)
saveRDS(evaluation, "artifacts/multiscale_emulator_evaluation.rds")
print(evaluation)
