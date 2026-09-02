#' Sample synthetic biological parameters for one host profile
#'
#' @param site_id Synthetic site identifier.
#' @param profile_id Profile identifier.
#' @param config Multiscale configuration.
#' @param seed Optional random seed.
#'
#' @return A one-row data frame containing profile metadata and truth parameters.
#' @export
sample_profile_parameters <- function(
    site_id,
    profile_id,
    config,
    seed = NULL
) {
  sites <- multiscale_site_table()
  site <- sites[match(site_id, sites$site_id), , drop = FALSE]
  if (!nrow(site)) stop("Unknown site_id: ", site_id)
  if (!is.null(seed)) set.seed(seed)
  multipliers <- c(
    h = site$h_multiplier,
    gamma = site$gamma_multiplier,
    delta = site$delta_multiplier
  )
  sampled <- exp(stats::rnorm(
    3L,
    log(config$truth_center * multipliers),
    config$truth_log_sd
  ))
  names(sampled) <- c("h", "gamma", "delta")
  bounds <- rbind(
    h = config$h_bounds,
    gamma = config$gamma_bounds,
    delta = config$delta_bounds
  )
  sampled <- pmin(pmax(sampled, bounds[, 1]), bounds[, 2])
  data.frame(
    profile_id = profile_id,
    site_id = site_id,
    h = sampled[["h"]],
    gamma = sampled[["gamma"]],
    delta = sampled[["delta"]],
    antibiotic_end = config$antibiotic_end,
    contact_rate = site$contact_rate,
    susceptible_fraction = site$susceptible_fraction,
    genomic_available = site$genomic_available,
    quantitative_available = site$quantitative_available,
    stringsAsFactors = FALSE
  )
}

#' Simulate noisy omic observations from a within-host trajectory
#'
#' @param trajectory Within-host trajectory.
#' @param profile One-row profile metadata from [sample_profile_parameters()].
#' @param config Multiscale configuration.
#' @param seed Optional random seed.
#'
#' @return A data frame with one row per scheduled observation.
#' @export
simulate_omic_observations <- function(
    trajectory,
    profile,
    config,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  rows <- match(config$observation_days, trajectory$time)
  if (anyNA(rows)) stop("trajectory must contain every observation day.")
  states <- trajectory[rows, , drop = FALSE]
  backgrounds <- config$backgrounds$background
  s_cols <- paste0("S_", backgrounds)
  r_cols <- paste0("R_", backgrounds)
  resistant <- as.matrix(states[r_cols])
  total_resistant <- rowSums(resistant)
  qpcr_truth <- log10(pmax(total_resistant * config$qpcr_scale, 1e-12))
  qpcr_observed <- qpcr_truth + stats::rnorm(
    length(qpcr_truth),
    0,
    config$qpcr_sd
  )
  qpcr_censored <- qpcr_observed < config$qpcr_lod
  qpcr_observed[qpcr_censored] <- config$qpcr_lod

  all_states <- as.matrix(states[c(s_cols, r_cols)])
  ecoli_cols <- c(1L, 2L, 5L, 6L)
  kleb_cols <- c(3L, 4L, 7L, 8L)
  ecoli_truth <- log10(pmax(rowSums(all_states[, ecoli_cols, drop = FALSE]), 1e-12))
  kleb_truth <- log10(pmax(rowSums(all_states[, kleb_cols, drop = FALSE]), 1e-12))
  ecoli_observed <- ecoli_truth + stats::rnorm(
    length(ecoli_truth),
    0,
    config$abundance_sd
  )
  kleb_observed <- kleb_truth + stats::rnorm(
    length(kleb_truth),
    0,
    config$abundance_sd
  )

  result <- data.frame(
    profile_id = profile$profile_id,
    site_id = profile$site_id,
    day = config$observation_days,
    antibiotic_active = antibiotic_activity(
      config$observation_days,
      config$antibiotic_start,
      profile$antibiotic_end
    ),
    qpcr_log10 = qpcr_observed,
    qpcr_censored = qpcr_censored,
    ecoli_log10 = ecoli_observed,
    klebsiella_log10 = kleb_observed,
    culture_positive = total_resistant * config$qpcr_scale >=
      10^config$qpcr_lod,
    stringsAsFactors = FALSE
  )
  for (index in seq_along(backgrounds)) {
    abundance <- log10(pmax(resistant[, index] * config$qpcr_scale, 1e-12))
    probability <- stats::plogis(
      config$linkage_logit_intercept +
        config$linkage_logit_slope * (abundance - config$qpcr_lod)
    )
    probability[resistant[, index] <= 0] <- 1 - config$linkage_specificity
    result[[paste0("link_", backgrounds[index])]] <- stats::rbinom(
      nrow(result),
      1L,
      probability
    )
  }
  if (!profile$quantitative_available) {
    result[c("qpcr_log10", "ecoli_log10", "klebsiella_log10")] <- NA_real_
  }
  if (!profile$genomic_available) {
    result[paste0("link_", backgrounds)] <- NA_integer_
  }
  result
}

#' Simulate one truth profile and its observation layer
#'
#' @param site_id Synthetic site identifier.
#' @param profile_id Profile identifier.
#' @param config Multiscale configuration.
#' @param seed Optional random seed.
#'
#' @return A list containing profile metadata, truth trajectory, and observations.
#' @export
simulate_multiscale_profile <- function(
    site_id,
    profile_id,
    config,
    seed = NULL
) {
  if (is.null(seed)) {
    seed <- as.integer(config$base_seed + sum(utf8ToInt(profile_id)))
  }
  profile <- sample_profile_parameters(site_id, profile_id, config, seed)
  parameters <- unlist(profile[c("h", "gamma", "delta")], use.names = TRUE)
  trajectory <- simulate_within_host(
    parameters,
    config,
    antibiotic_end = profile$antibiotic_end
  )
  observations <- simulate_omic_observations(
    trajectory,
    profile,
    config,
    seed + 1L
  )
  list(profile = profile, truth = trajectory, observations = observations)
}

#' Objective for fitting one within-host profile
#'
#' @param parameters Numeric vector ordered as `h`, `gamma`, and `delta`.
#' @param observations Longitudinal omic observations.
#' @param config Multiscale configuration.
#' @param antibiotic_end Treatment end time.
#'
#' @return Scalar weighted objective value.
#' @export
within_host_fit_objective <- function(
    parameters,
    observations,
    config,
    antibiotic_end = config$antibiotic_end
) {
  names(parameters) <- c("h", "gamma", "delta")
  trajectory <- tryCatch(
    simulate_within_host(
      parameters,
      config,
      antibiotic_end = antibiotic_end,
      end_time = max(observations$day),
      output_times = observations$day,
      step = config$fit_step
    ),
    error = function(error) NULL
  )
  if (is.null(trajectory)) return(1e12)
  backgrounds <- config$backgrounds$background
  susceptible <- as.matrix(trajectory[paste0("S_", backgrounds)])
  resistant <- as.matrix(trajectory[paste0("R_", backgrounds)])
  predicted_qpcr <- log10(pmax(rowSums(resistant) * config$qpcr_scale, 1e-12))
  ecoli <- rowSums(susceptible[, 1:2, drop = FALSE] + resistant[, 1:2, drop = FALSE])
  klebsiella <- rowSums(susceptible[, 3:4, drop = FALSE] + resistant[, 3:4, drop = FALSE])
  objective <- 0
  quantitative <- list(
    qpcr_log10 = cbind(predicted_qpcr, config$qpcr_sd),
    ecoli_log10 = cbind(log10(pmax(ecoli, 1e-12)), config$abundance_sd),
    klebsiella_log10 = cbind(log10(pmax(klebsiella, 1e-12)), config$abundance_sd)
  )
  for (field in names(quantitative)) {
    observed <- observations[[field]]
    keep <- is.finite(observed)
    if (any(keep)) {
      predicted <- quantitative[[field]][, 1]
      scale <- quantitative[[field]][, 2]
      objective <- objective + sum(((observed[keep] - predicted[keep]) / scale[keep])^2)
    }
  }
  for (index in seq_along(backgrounds)) {
    field <- paste0("link_", backgrounds[index])
    observed <- observations[[field]]
    keep <- !is.na(observed)
    if (any(keep)) {
      abundance <- log10(pmax(resistant[, index] * config$qpcr_scale, 1e-12))
      probability <- stats::plogis(
        config$linkage_logit_intercept +
          config$linkage_logit_slope * (abundance - config$qpcr_lod)
      )
      probability <- pmin(pmax(probability, 1e-8), 1 - 1e-8)
      objective <- objective - 2 * sum(stats::dbinom(
        observed[keep],
        1L,
        probability[keep],
        log = TRUE
      ))
    }
  }
  objective
}

#' Fit one within-host profile with bounded multi-start optimization
#'
#' @param observations Longitudinal observations for one profile.
#' @param config Multiscale configuration.
#' @param antibiotic_end Treatment end time.
#' @param n_starts Number of optimization starts.
#' @param seed Optional random seed.
#'
#' @return A list containing fitted parameters and optimization diagnostics.
#' @export
fit_within_host_profile <- function(
    observations,
    config,
    antibiotic_end = config$antibiotic_end,
    n_starts       = 5L,
    seed           = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  lower <- c(config$h_bounds[1], config$gamma_bounds[1], config$delta_bounds[1])
  upper <- c(config$h_bounds[2], config$gamma_bounds[2], config$delta_bounds[2])
  starts <- matrix(config$truth_center, nrow = 1L)
  if (n_starts > 1L) {
    random_starts <- t(replicate(
      n_starts - 1L,
      stats::runif(3L, lower, upper)
    ))
    starts <- rbind(starts, random_starts)
  }
  fits <- lapply(seq_len(n_starts), function(index) {
    stats::optim(
      starts[index, ],
      within_host_fit_objective,
      observations = observations,
      config = config,
      antibiotic_end = antibiotic_end,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 150L)
    )
  })
  objectives <- vapply(fits, `[[`, numeric(1), "value")
  best <- fits[[which.min(objectives)]]
  parameters <- stats::setNames(best$par, c("h", "gamma", "delta"))
  tolerance <- 1e-5
  list(
    parameters = parameters,
    objective = best$value,
    convergence = best$convergence,
    at_boundary = any(
      abs(parameters - lower) < tolerance |
        abs(parameters - upper) < tolerance
    ),
    starts = n_starts
  )
}
