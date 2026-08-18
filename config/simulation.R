default_simulation_config <- function() {
  list(
    n_agents = 1000L,
    prevalence = 0.01,
    contact_rate = 2,
    incubation_days = 3,
    recovery_rate = 0.20,
    beta_x = 1,
    alpha = -2.263589,
    delta_scenario = 0.352420,
    target_r0 = c(lower = 1.5, higher = 2.0),
    max_days = 250L,
    early_susceptible_fraction = 0.90,
    base_seed = 20260818L
  )
}

scenario_table <- function() {
  data.frame(
    scenario = c("lower", "higher"),
    scenario_high = c(0, 1),
    organism = c("organism_lower", "organism_higher"),
    stringsAsFactors = FALSE
  )
}

