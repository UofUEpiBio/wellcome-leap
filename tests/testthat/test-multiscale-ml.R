testthat::test_that("all seven multiscale modality patterns encode indicators", {
  testthat::skip_if_not_installed("torch")
  data <- data.frame(
    profile_id = "p1",
    site_id = "site_a",
    qpcr_baseline = 5,
    qpcr_peak = 6,
    qpcr_day30 = 4,
    ecoli_day30 = -0.2,
    klebsiella_day30 = -0.4,
    linked_backgrounds = 2,
    linkage_observations = 20,
    antibiotic_days = 7,
    contact_rate = 0.2,
    susceptible_fraction = 0.9,
    conjugation_multiplier = 1
  )
  preprocessor <- fit_multiscale_preprocessor(data)
  patterns <- multiscale_observation_patterns()
  testthat::expect_length(patterns, 7L)
  for (pattern in names(patterns)) {
    matrix <- make_multiscale_input_matrix(
      data,
      patterns[[pattern]],
      preprocessor
    )
    observed <- matrix[1, grep("_observed$", colnames(matrix))]
    testthat::expect_equal(sum(observed), length(patterns[[pattern]]))
  }
})

testthat::test_that("leave-site-out split prevents site leakage", {
  data <- data.frame(
    site_id = rep(paste0("site_", letters[1:4]), each = 3L),
    profile_id = seq_len(12L)
  )
  split <- split_multiscale_by_site(data)
  testthat::expect_setequal(unique(split$test$site_id), "site_d")
  testthat::expect_setequal(unique(split$validation$site_id), "site_c")
  testthat::expect_false(any(split$train$site_id %in% c("site_c", "site_d")))
})

testthat::test_that("profile holdout keeps scenarios together and every site represented", {
  data <- expand.grid(
    profile_number = seq_len(5L),
    site_id = paste0("site_", letters[1:4]),
    intervention = c("baseline", "shorter_antibiotic", "conjugation_inhibition"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$profile_id <- paste(data$site_id, data$profile_number)
  split <- split_multiscale_by_profile(data, seed = 92L)
  memberships <- lapply(split, function(partition) unique(partition$profile_id))
  testthat::expect_length(intersect(memberships$train, memberships$validation), 0L)
  testthat::expect_length(intersect(memberships$train, memberships$test), 0L)
  testthat::expect_length(intersect(memberships$validation, memberships$test), 0L)
  for (partition in split) {
    testthat::expect_setequal(unique(partition$site_id), paste0("site_", letters[1:4]))
    testthat::expect_true(all(table(partition$profile_id) == 3L))
  }
})

testthat::test_that("multiscale emulator trains and predicts positive metrics", {
  testthat::skip_if_not_installed("torch")
  set.seed(91L)
  data <- expand.grid(
    site_id = paste0("site_", letters[1:4]),
    profile_number = seq_len(4L),
    intervention = c("baseline", "shorter_antibiotic", "conjugation_inhibition"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$profile_id <- paste(data$site_id, data$profile_number)
  feature_names <- c(
    "qpcr_baseline", "qpcr_peak", "qpcr_day30", "ecoli_day30",
    "klebsiella_day30", "linked_backgrounds", "linkage_observations",
    "antibiotic_days", "contact_rate", "susceptible_fraction",
    "conjugation_multiplier"
  )
  for (field in feature_names) data[[field]] <- stats::runif(nrow(data))
  for (target in c("r0_within", "re_within", "r0_between", "re_between")) {
    data[[target]] <- exp(stats::rnorm(nrow(data), 0, 0.2))
  }
  fit <- fit_multiscale_emulator(
    data,
    hidden_dim_1 = 8L,
    hidden_dim_2 = 4L,
    max_epochs = 10L,
    patience = 3L,
    seed = 91L
  )
  prediction <- predict_multiscale_emulator(fit, fit$split$test)
  testthat::expect_true(all(as.matrix(prediction) > 0))
  testthat::expect_equal(nrow(prediction), nrow(fit$split$test))

  weights <- tempfile(fileext = ".pt")
  metadata <- tempfile(fileext = ".rds")
  save_multiscale_emulator(fit, weights, metadata)
  loaded <- load_multiscale_emulator(weights, metadata)
  loaded_prediction <- predict_multiscale_emulator(loaded, fit$split$test)
  testthat::expect_equal(prediction, loaded_prediction, tolerance = 1e-6)
})
