testthat::test_that("the web export reproduces torch predictions without torch", {
  testthat::skip_if_not_installed("torch")
  torch::torch_manual_seed(11L)
  training <- data.frame(
    x = seq(-2, 2, length.out = 40),
    scenario = rep(c("lower", "higher"), length.out = 40)
  )
  preprocessor <- fit_preprocessor(training)
  fit <- list(
    model = create_masked_count_model(hidden_dim_1 = 6L, hidden_dim_2 = 4L),
    metadata = list(
      preprocessor = preprocessor,
      architecture = list(input_dim = 4L, hidden_dim_1 = 6L, hidden_dim_2 = 4L),
      supported_patterns = c("both", "x_only", "scenario_only")
    )
  )
  web_model <- extract_web_model(fit)

  for (pattern in fit$metadata$supported_patterns) {
    input <- make_masked_matrix(training, pattern, preprocessor)
    testthat::expect_equal(
      predict_web_model(web_model, input),
      predict_rate_matrix(fit$model, input),
      tolerance = 1e-5
    )
  }
})

testthat::test_that("the web export writes the shapes the browser expects", {
  testthat::skip_if_not_installed("torch")
  testthat::skip_if_not_installed("jsonlite")
  training <- data.frame(x = c(-1, 1), scenario = c("lower", "higher"))
  fit <- list(
    model = create_masked_count_model(hidden_dim_1 = 6L, hidden_dim_2 = 4L),
    metadata = list(
      preprocessor = fit_preprocessor(training),
      architecture = list(input_dim = 4L, hidden_dim_1 = 6L, hidden_dim_2 = 4L),
      supported_patterns = c("both", "x_only", "scenario_only")
    )
  )
  path <- file.path(tempdir(), "web_model.json")
  on.exit(unlink(path), add = TRUE)
  export_web_model(fit, path)
  exported <- jsonlite::fromJSON(path, simplifyVector = FALSE)

  units <- c(hidden_1 = 6L, hidden_2 = 4L, output = 1L)
  inputs <- c(hidden_1 = 4L, hidden_2 = 6L, output = 4L)
  for (layer in names(units)) {
    testthat::expect_length(exported$layers[[layer]]$bias, units[[layer]])
    testthat::expect_length(
      exported$layers[[layer]]$weight,
      units[[layer]] * inputs[[layer]]
    )
  }
  testthat::expect_named(
    exported$preprocessor$scenario_levels,
    c("lower", "higher")
  )
})

testthat::test_that("the shipped export matches the shipped scenarios", {
  testthat::skip_if_not_installed("jsonlite")
  model_path <- file.path("..", "..", "app", "model.json")
  site_path <- file.path("..", "..", "app", "site.json")
  testthat::skip_if_not(file.exists(model_path) && file.exists(site_path))
  model <- jsonlite::fromJSON(model_path, simplifyVector = FALSE)
  site <- jsonlite::fromJSON(site_path, simplifyVector = FALSE)
  offered <- vapply(site$scenarios, function(scenario) scenario$id, character(1))
  testthat::expect_true(
    all(offered %in% names(model$preprocessor$scenario_levels))
  )
})
