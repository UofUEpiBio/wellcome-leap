#' Default simulation configuration
#'
#' @return A named list of epidemiological, calibration, and reproducibility
#'   settings for the toy experiment.
#' @export
default_simulation_config <- function() {
  list(
    n_agents                  = 10000L,
    prevalence                = 5 / 10000,
    contact_rate              = 4,
    incubation_days           = 7,
    recovery_rate             = 0.20,
    x_mean                    = 0,
    x_sd                      = 1,
    beta_x                    = 1,
    alpha                     = -4.124644,
    delta_scenario            = 2.474928,
    target_r0                 = c(lower = 0.5, higher = 4.0),
    max_days                  = 60L,
    early_susceptible_fraction = 0.90,
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
