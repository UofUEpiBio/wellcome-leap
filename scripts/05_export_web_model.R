# Export the fitted masked surrogate for the static browser application.
#
# Reads the generated torch weights and writes app/model.json, the only
# committed representation of the fitted model. Run this after
# scripts/04_fit_and_evaluate.R whenever the surrogate is retrained.

source("R/masking.R")
source("R/torch_model.R")
source("R/fit_models.R")
source("R/export_web_model.R")

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
