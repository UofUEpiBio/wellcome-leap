# Export the fitted masked surrogate for the static browser application.
#
# Reads the generated torch weights and writes app/model.json, the only
# committed representation of the fitted model. Run this after
# scripts/04_fit_and_evaluate.R whenever the surrogate is retrained.

source("R/masking.R")
source("R/torch_model.R")
source("R/fit_models.R")
source("R/export_web_model.R")
source("R/multimodal_emulator.R")

weights_path <- "artifacts/masked_model.pt"
metadata_path <- "artifacts/masked_model_metadata.rds"
rmse_path <- "artifacts/masked_model_rmse.rds"

if (!file.exists(weights_path) || !file.exists(metadata_path)) {
  stop("Run scripts/04_fit_and_evaluate.R before exporting the web model.")
}

fit <- load_masked_model(weights_path, metadata_path)
evaluation <- if (file.exists(rmse_path)) readRDS(rmse_path) else NULL
path <- export_web_model(fit, "app/model.json", evaluation = evaluation)

# Confirm the dependency-free forward pass reproduces torch before shipping.
web_model <- extract_web_model(fit)
grid <- expand.grid(
  x = seq(-3, 3, by = 0.25),
  scenario = c("lower", "higher"),
  stringsAsFactors = FALSE
)
patterns <- c("both", "x_only", "scenario_only")
largest_gap <- max(vapply(patterns, function(pattern) {
  input <- make_masked_matrix(grid, pattern, fit$metadata$preprocessor)
  max(abs(predict_rate_matrix(fit$model, input) -
            predict_web_model(web_model, input)))
}, numeric(1)))

cat("Wrote", path, "\n")
cat("Largest torch-versus-base-R gap:", format(largest_gap, digits = 3), "\n")
if (largest_gap > 1e-5) stop("Base R forward pass does not match torch.")

multiscale_weights <- "artifacts/multiscale_emulator.pt"
multiscale_metadata <- "artifacts/multiscale_emulator_metadata.rds"
multiscale_evaluation <- "artifacts/multiscale_emulator_evaluation.rds"
if (!file.exists(multiscale_weights) || !file.exists(multiscale_metadata)) {
  stop("Run scripts/08_train_multiscale_emulator.R before exporting the app.")
}
multiscale_fit <- load_multiscale_emulator(
  multiscale_weights,
  multiscale_metadata
)
multiscale_results <- if (file.exists(multiscale_evaluation)) {
  readRDS(multiscale_evaluation)
} else {
  NULL
}
export_multiscale_web_model(
  multiscale_fit,
  path,
  evaluation = multiscale_results
)
multiscale_web <- extract_multiscale_web_model(multiscale_fit)
multiscale_input <- matrix(
  seq(-1, 1, length.out = multiscale_fit$architecture$input_dim * 3L),
  nrow = 3L
)
multiscale_fit$model$eval()
torch::with_no_grad({
  multiscale_torch <- exp(as.matrix(multiscale_fit$model(torch::torch_tensor(
    multiscale_input,
    dtype = torch::torch_float()
  ))))
})
multiscale_gap <- max(abs(
  multiscale_torch -
    predict_multiscale_web_model(multiscale_web, multiscale_input)
))
cat("Appended multiscale weights to", path, "\n")
cat(
  "Largest multiscale torch-versus-base-R gap:",
  format(multiscale_gap, digits = 3),
  "\n"
)
if (multiscale_gap > 1e-5) {
  stop("Base R multiscale forward pass does not match torch.")
}
