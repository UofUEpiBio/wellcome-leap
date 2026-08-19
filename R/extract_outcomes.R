#' Exclude the epiworld seeding pseudo-source
#'
#' `epiworldR::run_multiple()` may represent initial infections with a
#' pseudo-source whose `source` value is `-1`. That record is not an agent's
#' realized reproduction number.
#'
#' @param reproduction Reproductive-number data containing a `source` column.
#'
#' @return The reproductive-number data restricted to nonnegative agent IDs.
#' @export
exclude_seed_pseudo_source <- function(reproduction) {
  if (!"source" %in% names(reproduction)) {
    stop("reproduction must contain a source column.")
  }
  reproduction[
    !is.na(reproduction$source) & reproduction$source >= 0,
    ,
    drop = FALSE
  ]
}

#' Last day of the approximately fully susceptible early phase
#'
#' An individual reproduction number only measures that individual's own
#' transmissibility while the susceptible pool is still essentially intact.
#' Once susceptibles are depleted, realized offspring counts are constrained by
#' the epidemic's final size rather than by the agent's transmission
#' probability.
#'
#' @param history Total-history data containing `date`, `state`, and `counts`.
#' @param config Named simulation-configuration list.
#'
#' @return The last day on which at least `early_susceptible_fraction` of the
#'   population was susceptible, or `-1` when the threshold is never met.
#' @export
early_phase_last_day <- function(
    history,
    config
) {
  required <- c("date", "state", "counts")
  if (!all(required %in% names(history))) {
    stop("history must contain: ", paste(required, collapse = ", "))
  }
  susceptible <- history[history$state == "Susceptible", , drop = FALSE]
  early_days <- susceptible$date[
    susceptible$counts / config$n_agents >= config$early_susceptible_fraction
  ]
  if (length(early_days)) max(early_days) else -1L
}

#' Flag individual reproduction numbers that are eligible for analysis
#'
#' An agent's realized `R_i` is analyzed only when the agent was given a full
#' opportunity to transmit and was infected before susceptible depletion
#' distorts offspring counts. Concretely, the agent must have been infected on
#' or before `infection_cutoff_day`, must have had at least
#' `min_transmission_days` of simulated time left after infection, and must
#' have been infected during the early, approximately fully susceptible phase.
#'
#' @param infection_day Integer vector of agent infection (exposure) days.
#' @param config Named simulation-configuration list.
#' @param early_last_day Last early-phase day, as returned by
#'   [early_phase_last_day()]. Use `Inf` to skip the depletion criterion.
#'
#' @return A logical vector flagging analyzable individual reproduction numbers.
#' @export
individual_ri_eligible <- function(
    infection_day,
    config,
    early_last_day = Inf
) {
  transmission_window <- config$max_days - infection_day
  !is.na(infection_day) &
    infection_day <= config$infection_cutoff_day &
    transmission_window >= config$min_transmission_days &
    infection_day <= early_last_day
}

#' Restrict agent outcomes to analysis-eligible individual reproduction numbers
#'
#' Applies whichever of the `outcome_complete` and `analysis_eligible` flags are
#' present, so the same filter can be used on production agent tables and on the
#' small synthetic tables used by the tests.
#'
#' @param agents Agent-level outcome data.
#'
#' @return The subset of `agents` whose individual `R_i` is analyzable.
#' @export
filter_analysis_agents <- function(agents) {
  keep <- rep(TRUE, nrow(agents))
  for (flag in c("outcome_complete", "analysis_eligible")) {
    if (flag %in% names(agents)) keep <- keep & agents[[flag]]
  }
  agents[keep, , drop = FALSE]
}

#' Extract individual and run-level transmission outcomes
#'
#' @param model_bundle List returned by [build_seirconn_model()].
#' @param config Named simulation-configuration list.
#' @param run_id Unique simulation-run identifier.
#' @param seed Integer random seed used to construct the run.
#'
#' @return A list containing agent outcomes, a run summary, and transmission
#'   edges.
#' @export
extract_simulation_results <- function(
    model_bundle,
    config,
    run_id,
    seed
) {
  model <- model_bundle$model
  agent_data <- model_bundle$agent_data
  reproduction <- epiworldR::get_reproductive_number(model)
  reproduction <- exclude_seed_pseudo_source(reproduction)
  transmissions <- epiworldR::get_transmissions(model)
  history <- epiworldR::get_hist_total(model)
  final_states <- epiworldR::get_agents_states(model)

  nonseed_edges <- transmissions[transmissions$source >= 0, , drop = FALSE]
  if (any(reproduction$source < 0)) {
    stop("Pseudo-source seeding rows must not enter individual outcomes.")
  }
  if (sum(reproduction$rt) != nrow(nonseed_edges)) {
    stop("Reproductive-number counts do not match non-seed transmission edges.")
  }

  final_day <- max(history$date)
  final_history <- history[history$date == final_day, , drop = FALSE]
  active_end <- sum(
    final_history$counts[final_history$state %in% c("Exposed", "Infected")]
  )
  run_complete <- active_end == 0

  early_last_day <- early_phase_last_day(history, config)

  infection_edge <- transmissions[
    match(reproduction$source, transmissions$target),
    c("target", "source", "date"),
    drop = FALSE
  ]
  idx <- match(reproduction$source, agent_data$agent_id)
  n_infected <- nrow(reproduction)
  agent_final_state <- final_states[reproduction$source + 1L]
  agent_outcome_complete <- !agent_final_state %in% c("Exposed", "Infected")
  infection_day <- reproduction$source_exposure_date
  agents <- data.frame(
    run_id = rep(run_id, n_infected),
    seed = rep(seed, n_infected),
    scenario = agent_data$scenario[idx],
    organism = agent_data$organism[idx],
    agent_id = reproduction$source,
    x = agent_data$x[idx],
    p_transmit = agent_data$p_transmit[idx],
    infection_day = infection_day,
    source_id = infection_edge$source,
    secondary_cases = reproduction$rt,
    final_state = agent_final_state,
    outcome_complete = agent_outcome_complete,
    transmission_window = config$max_days - infection_day,
    early_phase = infection_day <= early_last_day,
    analysis_eligible = agent_outcome_complete &
      individual_ri_eligible(infection_day, config, early_last_day),
    final_epidemic_size = rep(n_infected, n_infected),
    extinct_early = rep(nrow(nonseed_edges) == 0L, n_infected),
    stringsAsFactors = FALSE
  )

  run <- data.frame(
    run_id = run_id,
    seed = seed,
    scenario = agent_data$scenario[1],
    organism = agent_data$organism[1],
    final_day = final_day,
    early_last_day = early_last_day,
    final_epidemic_size = nrow(reproduction),
    nonseed_transmissions = nrow(nonseed_edges),
    active_at_end = active_end,
    outcome_complete = run_complete,
    extinct_early = nrow(nonseed_edges) == 0L,
    eligible_agents = sum(agents$analysis_eligible),
    eligible_mean_ri = if (any(agents$analysis_eligible)) {
      mean(agents$secondary_cases[agents$analysis_eligible])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  list(agents = agents, runs = run, transmissions = transmissions)
}
