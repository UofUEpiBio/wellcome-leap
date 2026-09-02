#' Default multiscale prototype configuration
#'
#' Values are deliberately synthetic. Four genomic backgrounds are nested in
#' two species, while the three fitted biological parameters are shared across
#' backgrounds to keep the five-point observation design identifiable.
#'
#' @return A named list of within-host, assay, between-host, and reproducibility
#'   settings.
#' @export
default_multiscale_config <- function() {
  backgrounds <- data.frame(
    background = paste0("bg", seq_len(4L)),
    species = rep(c("ecoli", "klebsiella"), each = 2L),
    growth_rate = c(0.80, 0.76, 0.70, 0.66),
    fitness_cost = c(0.08, 0.10, 0.07, 0.09),
    initial_share = c(0.32, 0.28, 0.22, 0.18),
    stringsAsFactors = FALSE
  )
  omega <- outer(
    backgrounds$species,
    backgrounds$species,
    function(recipient, donor) ifelse(recipient == donor, 1, 0.35)
  )
  diag(omega) <- 0.75
  dimnames(omega) <- list(backgrounds$background, backgrounds$background)

  list(
    backgrounds = backgrounds,
    omega = omega,
    carrying_capacity = 1,
    initial_total = 0.92,
    initial_resistant = 0.025,
    initial_donor = "bg1",
    ode_start = 0,
    observation_end = 30,
    carriage_horizon = 120,
    ngm_horizon = 365,
    ode_step = 0.10,
    fit_step = 0.25,
    ngm_step = 0.25,
    observation_days = c(0, 3, 7, 14, 30),
    antibiotic_start = 0,
    antibiotic_end = 7,
    shorter_antibiotic_end = 3,
    h_bounds = c(0.03, 0.35),
    gamma_bounds = c(0.04, 0.16),
    delta_bounds = c(0.04, 0.25),
    truth_center = c(h = 0.13, gamma = 0.09, delta = 0.13),
    truth_log_sd = c(h = 0.28, gamma = 0.20, delta = 0.24),
    ngm_seed = 1e-6,
    within_establishment = c(bg1 = 0.70, bg2 = 0.65, bg3 = 0.60, bg4 = 0.55),
    between_establishment = 0.55,
    inoculum_kappa = 2.5,
    carriage_detection = 1e-5,
    qpcr_scale = 1e9,
    qpcr_sd = 0.18,
    abundance_sd = 0.10,
    qpcr_lod = 2.5,
    linkage_logit_intercept = -3.2,
    linkage_logit_slope = 1.2,
    linkage_specificity = 0.99,
    base_seed = 20260901L
  )
}

#' Synthetic-site configuration
#'
#' Sites differ in biological parameter distributions, contact patterns, and
#' structural modality availability so leave-site-out evaluation is meaningful.
#'
#' @return A four-row data frame of synthetic-site settings.
#' @export
multiscale_site_table <- function() {
  data.frame(
    site_id = paste0("site_", letters[1:4]),
    h_multiplier = c(0.80, 1.00, 1.20, 1.35),
    gamma_multiplier = c(1.15, 1.00, 0.90, 0.80),
    delta_multiplier = c(0.85, 1.10, 1.00, 1.25),
    contact_rate = c(0.18, 0.22, 0.28, 0.32),
    susceptible_fraction = c(0.96, 0.92, 0.88, 0.84),
    genomic_available = c(TRUE, TRUE, FALSE, TRUE),
    quantitative_available = c(TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}
