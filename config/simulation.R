#' Default simulation configuration
#'
#' `max_days` is the simulated horizon, `infection_cutoff_day` is the last
#' infection day whose individual reproduction number is analyzed, and
#' `min_transmission_days` is the shortest transmission window an analyzed
#' individual must have been given. Together they keep every analyzed `R_i`
#' free of right-censoring: an agent infected on `infection_cutoff_day` still
#' has `max_days - infection_cutoff_day` days left to infect others.
#'
#' @return A named list of epidemiological, calibration, and reproducibility
#'   settings for the toy experiment.
#' @export
default_simulation_config <- function() {
  list(
    n_agents                   = 10000L,
    prevalence                 = 5 / 10000,
    contact_rate               = 4,
    incubation_days            = 7,
    recovery_rate              = 0.20,
    x_mean                     = 0,
    x_sd                       = 1,
    beta_x                     = 0.5,
    alpha                      = -2.616672,
    delta_scenario             = 0.984734,
    target_r0                  = c(lower = 1.5, higher = 3.0),
    max_days                   = 120L,
    infection_cutoff_day       = 60L,
    min_transmission_days      = 30L,
    early_susceptible_fraction = 0.95,
    base_seed                  = 20260818L
  )
}

#' Scenario metadata
#'
#' @return A two-row data frame mapping scenario labels to numeric indicators
#'   and organism labels.
#' @export
scenario_table <- function() {
  data.frame(
    scenario      = c("lower", "higher"),
    scenario_high = c(0, 1),
    organism      = c("organism_lower", "organism_higher"),
    stringsAsFactors = FALSE
  )
}
