testthat::test_that("run-level splitting prevents leakage", {
  data <- expand.grid(
    run_id = sprintf("%s_%02d", rep(c("lower", "higher"), each = 6), 1:6),
    agent = 1:3,
    KEEP.OUT.ATTRS = FALSE
  )
  data$scenario <- sub("_.*", "", data$run_id)
  split <- split_by_run(data, seed = 42L)
  testthat::expect_length(intersect(split$run_ids$train, split$run_ids$test), 0L)
  testthat::expect_length(
    intersect(split$run_ids$validation, split$run_ids$test),
    0L
  )
})

testthat::test_that("mask matrices encode all three patterns", {
  data <- data.frame(x = c(-0.8, 0.8), scenario = c("lower", "higher"))
  preprocessor <- fit_preprocessor(data)
  both <- make_masked_matrix(data, "both", preprocessor)
  x_only <- make_masked_matrix(data, "x_only", preprocessor)
  scenario_only <- make_masked_matrix(data, "scenario_only", preprocessor)
  testthat::expect_true(all(both[, 3:4] == 1))
  testthat::expect_true(all(x_only[, 2] == 0 & x_only[, 4] == 0))
  testthat::expect_true(all(scenario_only[, 1] == 0 & scenario_only[, 3] == 0))
})

testthat::test_that("scenario balancing equalizes the two scenarios", {
  scenario <- rep(c("lower", "higher"), times = c(10L, 90L))
  weights <- scenario_balance_weights(scenario)
  totals <- tapply(weights, scenario, sum)

  testthat::expect_equal(totals[["lower"]], totals[["higher"]])
  testthat::expect_equal(mean(weights), 1)
  testthat::expect_true(all(weights > 0))
})

testthat::test_that("torch modality dropout retains one or both inputs", {
  testthat::skip_if_not_installed("torch")
  torch::torch_manual_seed(19L)
  input <- torch::torch_tensor(
    matrix(rep(c(0.5, 1, 1, 1), 300L), ncol = 4L, byrow = TRUE),
    dtype = torch::torch_float()
  )
  dropped <- as.matrix(apply_modality_dropout(input))
  indicators <- unique(dropped[, 3:4, drop = FALSE])

  testthat::expect_true(all(rowSums(indicators) >= 1))
  testthat::expect_true(all(apply(indicators, 1, paste, collapse = "") %in%
    c("11", "10", "01")))
  testthat::expect_true(all(dropped[dropped[, 3] == 0, 1] == 0))
  testthat::expect_true(all(dropped[dropped[, 4] == 0, 2] == 0))
})

testthat::test_that("torch model trains, predicts, and serializes", {
  testthat::skip_if_not_installed("torch")
  set.seed(81L)
  runs <- expand.grid(
    scenario = c("lower", "higher"),
    replicate = 1:12,
    agent = 1:10,
    KEEP.OUT.ATTRS = FALSE
  )
  runs$run_id <- sprintf("%s_%02d", runs$scenario, runs$replicate)
  runs$x <- stats::rnorm(nrow(runs))
  rate <- exp(-0.5 + runs$x + 0.35 * (runs$scenario == "higher"))
  runs$secondary_cases <- stats::rpois(nrow(runs), rate)
  runs$outcome_complete <- TRUE

  fit <- fit_masked_model(
    runs,
    hidden_dim_1 = 8L,
    hidden_dim_2 = 4L,
    batch_size = 256L,
    max_epochs = 20L,
    patience = 5L,
    seed = 81L
  )
  testthat::expect_true(fit$metadata$balance_scenarios)
  one_row <- data.frame(x = 0.5, scenario = "higher")
  input <- make_masked_matrix(one_row, "both", fit$metadata$preprocessor)
  prediction <- predict_rate_matrix(fit$model, input)
  testthat::expect_true(is.finite(prediction))
  testthat::expect_gte(prediction, 0)

  weights <- tempfile(fileext = ".pt")
  metadata <- tempfile(fileext = ".rds")
  save_masked_model(fit, weights, metadata)
  loaded <- load_masked_model(weights, metadata)
  loaded_prediction <- predict_rate_matrix(loaded$model, input)
  testthat::expect_equal(
    prediction,
    loaded_prediction,
    tolerance = 1e-6
  )
})
