testthat::test_that("prediction wrapper supports all observation patterns", {
  testthat::skip_if_not_installed("torch")
  if (!exists("predict_secondary_cases")) testthat::skip("Wrapper not available.")
  set.seed(18L)
  data <- expand.grid(
    scenario = c("lower", "higher"),
    replicate = 1:8,
    agent = 1:8,
    KEEP.OUT.ATTRS = FALSE
  )
  data$run_id <- sprintf("%s_%02d", data$scenario, data$replicate)
  data$x <- stats::rnorm(nrow(data))
  data$secondary_cases <- stats::rpois(
    nrow(data),
    exp(-0.5 + data$x + 0.35 * (data$scenario == "higher"))
  )
  data$outcome_complete <- TRUE
  fit <- fit_masked_model(
    data,
    hidden_dim_1 = 8L,
    hidden_dim_2 = 4L,
    batch_size = 256L,
    max_epochs = 10L,
    patience = 3L,
    seed = 18L
  )

  testthat::expect_equal(
    predict_secondary_cases(fit, x = 0.5)$observation_pattern,
    "x_only"
  )
  testthat::expect_equal(
    predict_secondary_cases(fit, scenario = "lower")$observation_pattern,
    "scenario_only"
  )
  vector_prediction <- predict_secondary_cases(
    fit,
    x = c(0.2, NA, 0.8),
    scenario = c("lower", "higher", NA)
  )
  testthat::expect_equal(
    vector_prediction$observation_pattern,
    c("both", "scenario_only", "x_only")
  )
  testthat::expect_error(predict_secondary_cases(fit), "Provide x")
})
