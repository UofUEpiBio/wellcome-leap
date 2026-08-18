#' Validate simulation configuration
#'
#' @param config Named simulation-configuration list.
#'
#' @return The configuration invisibly; throws an error for invalid settings.
validate_simulation_config <- function(config) {
  required <- c(
    "n_agents", "prevalence", "contact_rate", "incubation_days",
    "recovery_rate", "x_mean", "x_sd", "beta_x", "alpha",
    "delta_scenario", "max_days", "early_susceptible_fraction", "base_seed"
  )
  missing <- setdiff(required, names(config))
  if (length(missing)) {
    stop("Missing configuration fields: ", paste(missing, collapse = ", "))
  }
  if (config$n_agents < 2L) stop("n_agents must be at least 2.")
  if (config$prevalence <= 0 || config$prevalence >= 1) {
    stop("prevalence must be strictly between 0 and 1.")
  }
  if (!is.finite(config$x_mean) || !is.finite(config$x_sd) || config$x_sd <= 0) {
    stop("x_mean must be finite and x_sd must be positive.")
  }
  probability_fields <- c("recovery_rate", "early_susceptible_fraction")
  for (field in probability_fields) {
    if (config[[field]] <= 0 || config[[field]] >= 1) {
      stop(field, " must be strictly between 0 and 1.")
    }
  }
  invisible(config)
}

#' Look up scenario metadata
#'
#' @param scenario Character scenario label.
#'
#' @return A one-row scenario metadata data frame.
scenario_info <- function(scenario) {
  scenarios <- scenario_table()
  idx <- match(scenario, scenarios$scenario)
  if (is.na(idx)) {
    stop("scenario must be one of: ", paste(scenarios$scenario, collapse = ", "))
  }
  scenarios[idx, , drop = FALSE]
}

#' Calculate agent-specific transmission probabilities
#'
#' @param x Numeric vector of finite agent features.
#' @param scenario Character scenario label.
#' @param config Named simulation-configuration list.
#'
#' @return A numeric vector of per-contact transmission probabilities.
#' @export
transmission_probability <- function(
    x,
    scenario,
    config
) {
  info <- scenario_info(scenario)
  stats::plogis(
    config$alpha +
      config$beta_x * x +
      config$delta_scenario * info$scenario_high
  )
}

#' Build a heterogeneous connected SEIR model
#'
#' @param config Named simulation-configuration list.
#' @param scenario Character scenario label.
#' @param x Numeric vector containing one feature value per agent.
#'
#' @return A list containing an `epiworldR` model pointer and agent metadata.
#' @export
build_seirconn_model <- function(
    config,
    scenario,
    x
) {
  validate_simulation_config(config)
  info <- scenario_info(scenario)
  if (length(x) != config$n_agents) {
    stop("x must contain one value per agent.")
  }
  if (any(!is.finite(x))) {
    stop("x must contain only finite values.")
  }
  if (!requireNamespace("epiworldR", quietly = TRUE)) {
    stop("The epiworldR package is required. Run scripts/00_install_dependencies.R.")
  }

  agents <- cbind(
    intercept = 1,
    x = as.numeric(x),
    scenario_high = rep(info$scenario_high, config$n_agents)
  )

  model <- epiworldR::ModelSEIRCONN(
    name = info$organism,
    n = config$n_agents,
    prevalence = config$prevalence,
    contact_rate = config$contact_rate,
    transmission_rate = stats::plogis(config$alpha),
    incubation_days = config$incubation_days,
    recovery_rate = config$recovery_rate
  )
  epiworldR::verbose_off(model)
  epiworldR::set_agents_data(model, agents)

  transmission_fun <- epiworldR::virus_fun_logit(
    vars = c(0L, 1L, 2L),
    coefs = c(config$alpha, config$beta_x, config$delta_scenario),
    model = model
  )
  epiworldR::set_prob_infecting_fun(
    virus = epiworldR::get_virus(model, 0),
    model = model,
    vfun = transmission_fun
  )

  list(
    model = model,
    agent_data = data.frame(
      agent_id = seq_len(config$n_agents) - 1L,
      x = x,
      scenario = scenario,
      scenario_high = info$scenario_high,
      organism = info$organism,
      p_transmit = transmission_probability(x, scenario, config),
      stringsAsFactors = FALSE
    )
  )
}
