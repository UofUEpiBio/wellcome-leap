testthat::test_that("breakpoint-aware ODE integration is stable and nonnegative", {
  config <- default_multiscale_config()
  parameters <- config$truth_center
  coarse <- simulate_within_host(
    parameters,
    config,
    end_time = 14,
    output_times = c(0, 3, 7, 14),
    step = 0.2
  )
  fine <- simulate_within_host(
    parameters,
    config,
    end_time = 14,
    output_times = c(0, 3, 7, 14),
    step = 0.1
  )
  testthat::expect_true(all(as.matrix(fine[-1]) >= 0))
  testthat::expect_equal(coarse, fine, tolerance = 1e-3)
})

testthat::test_that("transfer and division-associated loss affect intended states", {
  config <- default_multiscale_config()
  no_transfer <- config$truth_center
  no_transfer[["h"]] <- 0
  trajectory <- simulate_within_host(
    no_transfer,
    config,
    end_time = 10,
    output_times = c(0, 10)
  )
  testthat::expect_equal(trajectory$R_bg2[2], 0)
  testthat::expect_equal(trajectory$R_bg3[2], 0)
  testthat::expect_equal(trajectory$R_bg4[2], 0)

  stopped <- config
  stopped$backgrounds$growth_rate[] <- 0
  state <- initial_within_host_state(stopped)
  low <- within_host_derivative(
    1,
    state,
    c(h = 0, gamma = 0.04, delta = 0),
    stopped,
    antibiotic_end = 0
  )
  high <- within_host_derivative(
    1,
    state,
    c(h = 0, gamma = 0.16, delta = 0),
    stopped,
    antibiotic_end = 0
  )
  testthat::expect_identical(low, high)
})

testthat::test_that("synthetic observations preserve truth and site missingness", {
  config <- default_multiscale_config()
  available <- simulate_multiscale_profile("site_a", "available", config, 41L)
  no_genomics <- simulate_multiscale_profile("site_c", "no-genomics", config, 42L)
  no_quantitative <- simulate_multiscale_profile(
    "site_d",
    "no-quantitative",
    config,
    43L
  )
  testthat::expect_equal(available$observations$day, config$observation_days)
  testthat::expect_true(all(is.finite(available$truth$R_bg1)))
  testthat::expect_true(all(is.na(no_genomics$observations$link_bg1)))
  testthat::expect_true(all(is.na(no_quantitative$observations$qpcr_log10)))
})

testthat::test_that("noiseless observations recover the three shared parameters", {
  config <- default_multiscale_config()
  config$ode_step <- config$fit_step
  config$qpcr_sd <- 0.02
  config$abundance_sd <- 0.02
  parameters <- config$truth_center
  trajectory <- simulate_within_host(
    parameters,
    config,
    end_time = config$observation_end,
    output_times = config$observation_days
  )
  backgrounds <- config$backgrounds$background
  susceptible <- as.matrix(trajectory[paste0("S_", backgrounds)])
  resistant <- as.matrix(trajectory[paste0("R_", backgrounds)])
  observations <- data.frame(
    day = trajectory$time,
    qpcr_log10 = log10(rowSums(resistant) * config$qpcr_scale),
    ecoli_log10 = log10(rowSums(
      susceptible[, 1:2, drop = FALSE] + resistant[, 1:2, drop = FALSE]
    )),
    klebsiella_log10 = log10(rowSums(
      susceptible[, 3:4, drop = FALSE] + resistant[, 3:4, drop = FALSE]
    ))
  )
  for (background in backgrounds) {
    observations[[paste0("link_", background)]] <- NA_integer_
  }
  fit <- fit_within_host_profile(
    observations,
    config,
    n_starts = 1L,
    seed = 44L
  )
  testthat::expect_equal(fit$parameters, parameters, tolerance = 1e-3)
  testthat::expect_lt(fit$objective, 1e-6)
})

testthat::test_that("within-host next-generation matrix is a per-seed first generation", {
  config <- default_multiscale_config()
  parameters <- config$truth_center
  small <- within_host_next_generation(
    parameters,
    config,
    seed_size = 1e-6,
    horizon = 365
  )
  large <- within_host_next_generation(
    parameters,
    config,
    seed_size = 1e-5,
    horizon = 365
  )
  longer <- within_host_next_generation(
    parameters,
    config,
    seed_size = 1e-6,
    horizon = 500
  )
  testthat::expect_equal(small$lambda, large$lambda, tolerance = 1e-8)
  testthat::expect_equal(small$reproduction, longer$reproduction, tolerance = 2e-4)
  testthat::expect_true(all(diag(small$matrix) == 0))
  testthat::expect_true(all(small$matrix >= 0 & small$matrix <= 1))

  disabled <- parameters
  disabled[["h"]] <- 0
  zero <- within_host_next_generation(disabled, config)
  testthat::expect_equal(zero$reproduction, 0)
})

testthat::test_that("between-host Re integrates the trajectory and exposes Jensen bias", {
  config <- default_multiscale_config()
  metrics <- between_host_reproduction_metrics(
    config$truth_center,
    contact_rate = 0.25,
    susceptible_fraction = 0.9,
    config = config
  )
  effective <- metrics$re
  testthat::expect_gt(effective[["reproduction"]], 0)
  testthat::expect_gte(
    effective[["product_comparator"]],
    effective[["reproduction"]]
  )

  lower_contact <- between_host_reproduction(
    metrics$current_trajectory,
    0.10,
    0.9,
    config
  )
  testthat::expect_gt(
    effective[["reproduction"]],
    lower_contact[["reproduction"]]
  )
})

testthat::test_that("paired interventions change only declared mechanisms", {
  config <- default_multiscale_config()
  simulated <- simulate_multiscale_profile("site_b", "paired", config, 51L)
  bundle <- build_multiscale_profile_bundle(simulated, config)
  targets <- bundle$targets
  baseline <- targets[targets$intervention == "baseline", ]
  shorter <- targets[targets$intervention == "shorter_antibiotic", ]
  inhibited <- targets[targets$intervention == "conjugation_inhibition", ]
  testthat::expect_equal(baseline$r0_within, shorter$r0_within)
  testthat::expect_lt(shorter$re_between, baseline$re_between)
  testthat::expect_lt(inhibited$re_within, baseline$re_within)
})

testthat::test_that("summary-coupled ABM extraction agrees with transmission edges", {
  profile_data <- data.frame(
    profile_id = "fixture",
    effective_contact_risk = 0.20,
    carriage_duration = 20,
    contact_rate = 0.30,
    susceptible_fraction = 0.95
  )
  result <- run_multiscale_abm(
    profile_data,
    n_agents = 300L,
    prevalence = 0.03,
    days = 40L,
    seed = 202L
  )
  testthat::expect_equal(
    sum(result$agents$realized_secondary_acquisitions),
    sum(result$transmissions$source >= 0)
  )
})
