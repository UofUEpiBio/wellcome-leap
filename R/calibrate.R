log1pexp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

mean_logit_uniform <- function(intercept, beta = 1) {
  if (beta == 0) return(stats::plogis(intercept))
  (log1pexp(intercept + beta) - log1pexp(intercept)) / beta
}

analytic_effective_intercept <- function(target_r0, config) {
  target_mean_probability <- target_r0 * config$recovery_rate / config$contact_rate
  if (target_mean_probability <= 0 || target_mean_probability >= 1) {
    stop("The requested target implies an invalid mean transmission probability.")
  }
  stats::uniroot(
    function(intercept) {
      mean_logit_uniform(intercept, config$beta_x) - target_mean_probability
    },
    interval = c(-20, 20)
  )$root
}

analytic_starting_parameters <- function(config) {
  low <- analytic_effective_intercept(config$target_r0[["lower"]], config)
  high <- analytic_effective_intercept(config$target_r0[["higher"]], config)
  c(alpha = low, delta_scenario = high - low)
}

set_effective_intercept <- function(config, scenario, intercept) {
  if (scenario == "lower") {
    config$alpha <- intercept
  } else if (scenario == "higher") {
    config$delta_scenario <- intercept - config$alpha
  } else {
    stop("Unknown scenario: ", scenario)
  }
  config
}

estimate_early_reproduction <- function(config, scenario, n_reps = 200L,
                                        workers = 1L) {
  result <- run_simulation_study(
    config, n_reps = n_reps, scenarios = scenario, workers = workers
  )
  eligible <- result$agents$outcome_complete & result$agents$early_phase
  values <- result$agents$secondary_cases[eligible]
  c(
    mean = mean(values),
    se = stats::sd(values) / sqrt(length(values)),
    n_agents = length(values),
    complete_runs = sum(result$runs$outcome_complete)
  )
}

calibrate_scenario <- function(config, scenario, target_r0, n_reps = 200L,
                               iterations = 5L, workers = 1L,
                               half_width = 0.6) {
  center <- analytic_effective_intercept(target_r0, config)
  lower <- center - half_width
  upper <- center + half_width
  history <- vector("list", iterations)

  for (iteration in seq_len(iterations)) {
    candidate <- (lower + upper) / 2
    candidate_config <- set_effective_intercept(config, scenario, candidate)
    estimate <- estimate_early_reproduction(
      candidate_config, scenario, n_reps = n_reps, workers = workers
    )
    history[[iteration]] <- data.frame(
      iteration = iteration,
      scenario = scenario,
      target_r0 = target_r0,
      effective_intercept = candidate,
      estimated_r0 = estimate[["mean"]],
      se = estimate[["se"]],
      n_agents = estimate[["n_agents"]],
      stringsAsFactors = FALSE
    )
    if (estimate[["mean"]] < target_r0) lower <- candidate else upper <- candidate
  }

  history <- do.call(rbind, history)
  best <- which.min(abs(history$estimated_r0 - target_r0))
  list(
    effective_intercept = history$effective_intercept[best],
    history = history
  )
}

calibrate_scenarios <- function(config, n_reps = 200L, iterations = 5L,
                                workers = 1L) {
  low <- calibrate_scenario(
    config, "lower", config$target_r0[["lower"]], n_reps, iterations, workers
  )
  config$alpha <- low$effective_intercept
  high <- calibrate_scenario(
    config, "higher", config$target_r0[["higher"]], n_reps, iterations, workers
  )
  config$delta_scenario <- high$effective_intercept - config$alpha
  list(
    config = config,
    history = rbind(low$history, high$history)
  )
}

