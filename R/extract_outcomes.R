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

  susceptible <- history[history$state == "Susceptible", , drop = FALSE]
  early_days <- susceptible$date[
    susceptible$counts / config$n_agents >= config$early_susceptible_fraction
  ]
  early_last_day <- if (length(early_days)) max(early_days) else -1L

  infection_edge <- transmissions[
    match(reproduction$source, transmissions$target),
    c("target", "source", "date"),
    drop = FALSE
  ]
  idx <- match(reproduction$source, agent_data$agent_id)
  n_infected <- nrow(reproduction)
  agent_final_state <- final_states[reproduction$source + 1L]
  agent_outcome_complete <- !agent_final_state %in% c("Exposed", "Infected")
  agents <- data.frame(
    run_id = rep(run_id, n_infected),
    seed = rep(seed, n_infected),
    scenario = agent_data$scenario[idx],
    organism = agent_data$organism[idx],
    agent_id = reproduction$source,
    x = agent_data$x[idx],
    p_transmit = agent_data$p_transmit[idx],
    infection_day = reproduction$source_exposure_date,
    source_id = infection_edge$source,
    secondary_cases = reproduction$rt,
    final_state = agent_final_state,
    outcome_complete = agent_outcome_complete,
    early_phase = reproduction$source_exposure_date <= early_last_day,
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
    final_epidemic_size = nrow(reproduction),
    nonseed_transmissions = nrow(nonseed_edges),
    active_at_end = active_end,
    outcome_complete = run_complete,
    extinct_early = nrow(nonseed_edges) == 0L,
    early_mean_ri = if (any(agents$early_phase & agents$outcome_complete)) {
      mean(agents$secondary_cases[
        agents$early_phase & agents$outcome_complete
      ])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  list(agents = agents, runs = run, transmissions = transmissions)
}
