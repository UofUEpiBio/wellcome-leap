testthat::test_that("analytic starting parameters match requested targets", {
  config <- default_simulation_config()
  pars <- analytic_starting_parameters(config)
  testthat::expect_equal(unname(pars[["alpha"]]), -2.2636, tolerance = 1e-3)
  testthat::expect_gt(pars[["delta_scenario"]], 0)
})

testthat::test_that("transmission probability increases with x and scenario", {
  config <- default_simulation_config()
  x <- c(0, 0.5, 1)
  low <- transmission_probability(x, "lower", config)
  high <- transmission_probability(x, "higher", config)
  testthat::expect_true(all(diff(low) > 0))
  testthat::expect_true(all(high > low))
  testthat::expect_true(all(c(low, high) > 0 & c(low, high) < 1))
})

testthat::test_that("ModelSEIRCONN callback and extraction agree with edges", {
  config <- default_simulation_config()
  config$n_agents <- 200L
  config$max_days <- 100L
  result <- simulate_one(
    config, "lower", replicate_id = 1L, seed = 90210L,
    keep_transmissions = TRUE
  )
  testthat::expect_equal(
    sum(result$agents$secondary_cases),
    sum(result$transmissions$source >= 0)
  )
  testthat::expect_true(any(result$agents$secondary_cases == 0))
  testthat::expect_true(all(result$agents$x >= 0 & result$agents$x <= 1))
})

testthat::test_that("simulation is reproducible", {
  config <- default_simulation_config()
  config$n_agents <- 100L
  config$max_days <- 50L
  a <- simulate_one(config, "higher", 1L, seed = 123L)
  b <- simulate_one(config, "higher", 1L, seed = 123L)
  testthat::expect_identical(a$agents, b$agents)
  testthat::expect_identical(a$runs, b$runs)
})

