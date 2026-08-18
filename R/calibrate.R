#' Mean logistic probability over a normal feature
#'
#' @param intercept Logistic-model intercept.
#' @param beta Coefficient for the normal feature.
#' @param mean Feature-distribution mean.
#' @param sd Feature-distribution standard deviation.
#'
#' @return Expected logistic probability.
mean_logit_normal <- function(
    intercept,
    beta,
    mean = 0,
    sd   = 1
) {
  stats::integrate(
    function(x) {
      stats::plogis(intercept + beta * x) * stats::dnorm(x, mean, sd)
    },
    lower = -Inf,
    upper = Inf,
    subdivisions = 200L,
    rel.tol = 1e-9
  )$value
}

#' Solve an analytic starting intercept for a target reproduction number
#'
#' @param target_r0 Target early-epidemic reproduction number.
#' @param config Named simulation-configuration list.
#'
#' @return Numeric effective intercept.
#' @export
analytic_effective_intercept <- function(
    target_r0,
    config
) {
  target_mean_probability <- target_r0 * config$recovery_rate / config$contact_rate
  if (target_mean_probability <= 0 || target_mean_probability >= 1) {
    stop("The requested target implies an invalid mean transmission probability.")
  }
  stats::uniroot(
    function(intercept) {
      mean_logit_normal(
        intercept,
        beta = config$beta_x,
        mean = config$x_mean,
        sd = config$x_sd
      ) - target_mean_probability
    },
    interval = c(-20, 20)
  )$root
}

#' Calculate analytic scenario starting parameters
#'
#' @param config Named simulation-configuration list.
#'
#' @return Named numeric vector containing `alpha` and `delta_scenario`.
#' @export
analytic_starting_parameters <- function(config) {
  low <- analytic_effective_intercept(config$target_r0[["lower"]], config)
  high <- analytic_effective_intercept(config$target_r0[["higher"]], config)
  c(alpha = low, delta_scenario = high - low)
}

#' Set a scenario's effective transmission intercept
#'
#' @param config Named simulation-configuration list.
#' @param scenario Character scenario label.
#' @param intercept Numeric effective intercept.
#'
#' @return Updated configuration list.
set_effective_intercept <- function(
    config,
    scenario,
    intercept
) {
  if (scenario == "lower") {
    config$alpha <- intercept
  } else if (scenario == "higher") {
    config$delta_scenario <- intercept - config$alpha
  } else {
    stop("Unknown scenario: ", scenario)
  }
  config
}

#' Estimate early-epidemic reproduction in simulation
#'
#' @param config Named simulation-configuration list.
#' @param scenario Character scenario label.
#' @param n_reps Number of simulation replicates.
#' @param workers Number of parallel workers.
#'
#' @return Named numeric vector with the mean, standard error, agent count, and
#'   complete-run count.
#' @export
estimate_early_reproduction <- function(
    config,
    scenario,
    n_reps  = 200L,
    workers = 1L
) {
  result <- run_simulation_study(
    config,
    n_reps   = n_reps,
    scenarios = scenario,
    workers   = workers
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

#' Calibrate one scenario's transmission intercept
#'
#' @param config Named simulation-configuration list.
#' @param scenario Character scenario label.
#' @param target_r0 Target early-epidemic reproduction number.
#' @param n_reps Replicates per calibration candidate.
#' @param iterations Number of bounded-search iterations.
#' @param workers Number of parallel workers.
#' @param half_width Initial half-width around the analytic intercept.
#'
#' @return A list containing the selected effective intercept and search history.
#' @export
calibrate_scenario <- function(
    config,
    scenario,
    target_r0,
    n_reps     = 200L,
    iterations = 5L,
    workers    = 1L,
    half_width = 0.6
) {
  center <- analytic_effective_intercept(target_r0, config)
  lower <- center - half_width
  upper <- center + half_width
  history <- vector("list", iterations)

  for (iteration in seq_len(iterations)) {
    candidate <- (lower + upper) / 2
    candidate_config <- set_effective_intercept(config, scenario, candidate)
    estimate <- estimate_early_reproduction(
      candidate_config,
      scenario,
      n_reps = n_reps,
      workers = workers
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

#' Calibrate both organism scenarios
#'
#' @param config Named simulation-configuration list.
#' @param n_reps Replicates per calibration candidate.
#' @param iterations Number of bounded-search iterations.
#' @param workers Number of parallel workers.
#'
#' @return A list containing the calibrated configuration and combined history.
#' @export
calibrate_scenarios <- function(
    config,
    n_reps     = 200L,
    iterations = 5L,
    workers    = 1L
) {
  low <- calibrate_scenario(
    config,
    "lower",
    config$target_r0[["lower"]],
    n_reps,
    iterations,
    workers
  )
  config$alpha <- low$effective_intercept
  high <- calibrate_scenario(
    config,
    "higher",
    config$target_r0[["higher"]],
    n_reps,
    iterations,
    workers
  )
  config$delta_scenario <- high$effective_intercept - config$alpha
  list(
    config = config,
    history = rbind(low$history, high$history)
  )
}
