source("R/multimodal_emulator.R")
source("R/multiscale.R")

study_path <- "data/derived/multiscale/fitted_study.rds"
if (!file.exists(study_path)) {
  stop("Run scripts/07_fit_multiscale.R before training the emulator.")
}
study <- readRDS(study_path)
data <- make_multiscale_emulator_data(study$fitted)
site_fit <- fit_multiscale_emulator(
  data,
  hidden_dim_1 = 64L,
  hidden_dim_2 = 32L,
  learning_rate = 0.003,
  max_epochs = 500L,
  patience = 50L
)
site_evaluation <- evaluate_multiscale_emulator(site_fit)
site_evaluation$evaluation_scope <- "leave_site_out"

profile_partitions <- split_multiscale_by_profile(data)
fit <- fit_multiscale_emulator(
  data,
  hidden_dim_1 = 64L,
  hidden_dim_2 = 32L,
  learning_rate = 0.003,
  max_epochs = 500L,
  patience = 50L,
  partitions = profile_partitions
)
profile_evaluation <- evaluate_multiscale_emulator(fit)
profile_evaluation$evaluation_scope <- "profile_holdout"
evaluation <- rbind(site_evaluation, profile_evaluation)
evaluation <- evaluation[
  c(
    "evaluation_scope", "pattern", "target", "n", "rmse",
    "balanced_accuracy"
  )
]
save_multiscale_emulator(fit)
saveRDS(evaluation, "artifacts/multiscale_emulator_evaluation.rds")
print(evaluation)
