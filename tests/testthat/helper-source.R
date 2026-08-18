source(file.path("..", "..", "config", "simulation.R"))
source(file.path("..", "..", "R", "build_model.R"))
source(file.path("..", "..", "R", "extract_outcomes.R"))
source(file.path("..", "..", "R", "simulate.R"))
source(file.path("..", "..", "R", "calibrate.R"))
source(file.path("..", "..", "R", "masking.R"))
if (requireNamespace("torch", quietly = TRUE)) {
  source(file.path("..", "..", "R", "torch_model.R"))
  source(file.path("..", "..", "R", "fit_models.R"))
  predict_path <- file.path("..", "..", "R", "predict.R")
  if (file.exists(predict_path)) source(predict_path)
}
