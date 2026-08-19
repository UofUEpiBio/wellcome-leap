testthat::test_that("analytic starting parameters match requested targets", {
  config <- default_simulation_config()
  pars <- analytic_starting_parameters(config)
  effective_intercepts <- c(
    lower = pars[["alpha"]],
    higher = pars[["alpha"]] + pars[["delta_scenario"]]
  )
  implied_r0 <- vapply(
    effective_intercepts,
    mean_logit_normal,
    numeric(1),
    beta = config$beta_x,
    mean = config$x_mean,
    sd = config$x_sd
  ) * config$contact_rate / config$recovery_rate
  testthat::expect_equal(implied_r0, config$target_r0, tolerance = 1e-3)
  testthat::expect_gt(pars[["delta_scenario"]], 0)
})

testthat::test_that("transmission probability increases with x and scenario", {
  config <- default_simulation_config()
  x <- c(-1, 0, 1)
  low <- transmission_probability(x, "lower", config)
  high <- transmission_probability(x, "higher", config)
  testthat::expect_true(all(diff(low) > 0))
  testthat::expect_true(all(high > low))
  testthat::expect_true(all(c(low, high) > 0 & c(low, high) < 1))
})

testthat::test_that("seeding pseudo-source rows are excluded", {
  reproduction <- data.frame(
    source = c(-1L, 0L, 1L),
    rt = c(2L, 1L, 0L)
  )
  filtered <- exclude_seed_pseudo_source(reproduction)

  testthat::expect_equal(filtered$source, c(0L, 1L))
  testthat::expect_true(all(filtered$source >= 0L))
})

testthat::test_that("ModelSEIRCONN callback and extraction agree with edges", {
  config <- default_simulation_config()
  config$n_agents <- 200L
  config$prevalence <- 0.05
  config$max_days <- 100L
  result <- simulate_one(
    config,
    "higher",
    replicate_id      = 1L,
    seed              = 90210L,
    keep_transmissions = TRUE
  )
  testthat::expect_equal(
    sum(result$agents$secondary_cases),
    sum(result$transmissions$source >= 0)
  )
  testthat::expect_true(any(result$agents$secondary_cases == 0))
  testthat::expect_true(any(result$agents$x < 0))
  testthat::expect_true(any(result$agents$x > 0))
  testthat::expect_true(all(result$agents$agent_id >= 0))
})

testthat::test_that("simulation is reproducible", {
  config <- default_simulation_config()
  config$n_agents <- 100L
  config$max_days <- 50L
  config$infection_cutoff_day <- 20L
  a <- simulate_one(config, "higher", 1L, seed = 123L)
  b <- simulate_one(config, "higher", 1L, seed = 123L)
  testthat::expect_identical(a$agents, b$agents)
  testthat::expect_identical(a$runs, b$runs)
})

testthat::test_that("production batches can be reused and loaded", {
  config <- default_simulation_config()
  config$n_agents <- 50L
  config$max_days <- 40L
  config$infection_cutoff_day <- 10L
  output_dir <- tempfile("simulation-batches-")
  first_manifest <- run_simulation_batches(
    config,
    n_reps         = 3L,
    batch_size     = 2L,
    workers        = 1L,
    output_dir     = output_dir,
    reuse_existing = TRUE
  )
  second_manifest <- run_simulation_batches(
    config,
    n_reps         = 3L,
    batch_size     = 2L,
    workers        = 1L,
    output_dir     = output_dir,
    reuse_existing = TRUE
  )
  study <- load_simulation_batches(
    manifest_path = file.path(output_dir, "manifest.rds")
  )

  testthat::expect_identical(first_manifest, second_manifest)
  testthat::expect_equal(nrow(study$runs), 6L)
  testthat::expect_equal(length(unique(study$runs$run_id)), 6L)
})

testthat::test_that("the analysis window rejects censored configurations", {
  config <- default_simulation_config()
  testthat::expect_gte(
    config$max_days - config$infection_cutoff_day,
    config$min_transmission_days
  )

  censored <- config
  censored$infection_cutoff_day <- censored$max_days - 1L
  testthat::expect_error(
    validate_simulation_config(censored),
    "transmission window"
  )
})

testthat::test_that("individual Ri eligibility applies both cutoffs", {
  config <- default_simulation_config()
  infection_day <- c(0L, 60L, 61L)

  testthat::expect_equal(
    individual_ri_eligible(infection_day, config),
    c(TRUE, TRUE, FALSE)
  )
  testthat::expect_equal(
    individual_ri_eligible(infection_day, config, early_last_day = 30L),
    c(TRUE, FALSE, FALSE)
  )
  testthat::expect_true(all(
    config$max_days -
      infection_day[individual_ri_eligible(infection_day, config)] >=
      config$min_transmission_days
  ))
})

testthat::test_that("eligible agents exclude late and depleted infections", {
  config <- default_simulation_config()
  config$n_agents <- 500L
  config$prevalence <- 0.01
  result <- simulate_one(config, "higher", replicate_id = 3L, seed = 5150L)
  agents <- result$agents
  eligible <- filter_analysis_agents(agents)

  testthat::expect_true(all(
    eligible$infection_day <= config$infection_cutoff_day
  ))
  testthat::expect_true(all(
    eligible$transmission_window >= config$min_transmission_days
  ))
  testthat::expect_true(all(eligible$outcome_complete))
  testthat::expect_true(all(
    eligible$infection_day <= result$runs$early_last_day
  ))
  testthat::expect_equal(nrow(eligible), result$runs$eligible_agents)
  testthat::expect_equal(
    mean(eligible$secondary_cases),
    result$runs$eligible_mean_ri
  )
})
