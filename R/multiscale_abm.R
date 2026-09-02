#' Build a summary-coupled pARG carriage ABM
#'
#' @param profile_data Profile-level data containing effective contact risk,
#'   carriage duration, contact rate, and susceptible fraction.
#' @param n_agents Number of agents.
#' @param prevalence Initial carrier prevalence.
#' @param seed Random seed used to assign profiles.
#'
#' @return A list containing an `epiworldR` model and assigned agent profiles.
#' @export
build_multiscale_abm <- function(
    profile_data,
    n_agents  = 2000L,
    prevalence = 5 / 2000,
    seed       = 20260901L
) {
  if (!requireNamespace("epiworldR", quietly = TRUE)) {
    stop("The epiworldR package is required.")
  }
  required <- c(
    "profile_id", "effective_contact_risk", "carriage_duration",
    "contact_rate", "susceptible_fraction"
  )
  if (!all(required %in% names(profile_data))) {
    stop("profile_data is missing required summary columns.")
  }
  set.seed(seed)
  assigned <- profile_data[sample.int(
    nrow(profile_data),
    n_agents,
    replace = TRUE
  ), , drop = FALSE]
  base_contact <- max(assigned$contact_rate)
  probability <- assigned$effective_contact_risk *
    assigned$contact_rate / base_contact *
    assigned$susceptible_fraction * 0.55
  probability <- pmin(pmax(probability, 1e-8), 1 - 1e-8)
  recovery <- pmin(pmax(1 / pmax(assigned$carriage_duration, 1), 1e-8), 1 - 1e-8)
  agent_matrix <- cbind(
    transmission_logit = stats::qlogis(probability),
    recovery_logit = stats::qlogis(recovery)
  )
  model <- epiworldR::ModelSEIRCONN(
    name = "synthetic_pARG",
    n = n_agents,
    prevalence = prevalence,
    contact_rate = base_contact,
    transmission_rate = mean(probability),
    incubation_days = 1,
    recovery_rate = mean(recovery)
  )
  epiworldR::verbose_off(model)
  epiworldR::set_agents_data(model, agent_matrix)
  transmission_fun <- epiworldR::virus_fun_logit(
    vars = 0L,
    coefs = 1,
    model = model
  )
  recovery_fun <- epiworldR::virus_fun_logit(
    vars = 1L,
    coefs = 1,
    model = model
  )
  virus <- epiworldR::get_virus(model, 0)
  epiworldR::set_prob_infecting_fun(virus, model, transmission_fun)
  epiworldR::set_prob_recovery_fun(virus, model, recovery_fun)
  assigned$agent_id <- seq_len(n_agents) - 1L
  assigned$transmission_probability <- probability
  assigned$recovery_probability <- recovery
  list(model = model, agents = assigned)
}

#' Run a summary-coupled pARG carriage ABM
#'
#' @param profile_data Profile summaries accepted by [build_multiscale_abm()].
#' @param n_agents Number of agents.
#' @param prevalence Initial carrier prevalence.
#' @param days Simulation horizon.
#' @param seed Random seed.
#'
#' @return A list containing individual realized acquisitions, transmissions,
#'   and the assigned profile table.
#' @export
run_multiscale_abm <- function(
    profile_data,
    n_agents  = 2000L,
    prevalence = 5 / 2000,
    days       = 120L,
    seed       = 20260901L
) {
  built <- build_multiscale_abm(profile_data, n_agents, prevalence, seed)
  epiworldR::run(built$model, ndays = days, seed = seed + 1L)
  reproduction <- epiworldR::get_reproductive_number(built$model)
  reproduction <- reproduction[
    !is.na(reproduction$source) & reproduction$source >= 0,
    ,
    drop = FALSE
  ]
  index <- match(reproduction$source, built$agents$agent_id)
  agents <- cbind(
    built$agents[index, , drop = FALSE],
    infection_day = reproduction$source_exposure_date,
    realized_secondary_acquisitions = reproduction$rt
  )
  list(
    agents = agents,
    transmissions = epiworldR::get_transmissions(built$model),
    assigned_profiles = built$agents
  )
}
