#' Construct a seed-normalized first-generation within-host matrix
#'
#' Recipient acquisitions are not added back as donors. Matrix entries are the
#' probability that one introduced donor lineage successfully establishes the
#' pARG in a previously negative background.
#'
#' @param parameters Named `h`, `gamma`, and `delta` parameters.
#' @param config Multiscale configuration.
#' @param antibiotic_end Treatment end time.
#' @param seed_size Introduced donor abundance used for numerical integration.
#' @param horizon Integration horizon.
#'
#' @return A list containing exposure matrix `lambda`, next-generation matrix
#'   `matrix`, and its dominant eigenvalue `reproduction`.
#' @export
within_host_next_generation <- function(
    parameters,
    config,
    antibiotic_end = config$antibiotic_end,
    seed_size      = config$ngm_seed,
    horizon        = config$ngm_horizon
) {
  backgrounds <- config$backgrounds$background
  p_arg_free <- initial_within_host_state(config, resistant = FALSE)
  susceptible_trajectory <- simulate_within_host(
    parameters,
    config,
    initial_state = p_arg_free,
    antibiotic_end = antibiotic_end,
    end_time = horizon,
    output_times = sort(unique(c(
      seq(config$ode_start, horizon, by = config$ngm_step),
      horizon
    ))),
    step = config$ngm_step
  )
  time <- susceptible_trajectory$time
  s_matrix <- as.matrix(susceptible_trajectory[paste0("S_", backgrounds)])
  total_reference <- rowSums(s_matrix)
  lambda <- matrix(
    0,
    nrow = length(backgrounds),
    ncol = length(backgrounds),
    dimnames = list(recipient = backgrounds, donor = backgrounds)
  )
  for (donor in seq_along(backgrounds)) {
    growth <- config$backgrounds$growth_rate[donor] *
      (1 - config$backgrounds$fitness_cost[donor])
    donor_lineage <- numeric(length(time))
    donor_lineage[1] <- seed_size
    for (index in seq_len(length(time) - 1L)) {
      dt <- time[index + 1L] - time[index]
      rate_1 <- growth * (
        1 - parameters[["gamma"]] -
          total_reference[index] / config$carrying_capacity
      )
      rate_2 <- growth * (
        1 - parameters[["gamma"]] -
          mean(total_reference[index:(index + 1L)]) / config$carrying_capacity
      )
      rate_4 <- growth * (
        1 - parameters[["gamma"]] -
          total_reference[index + 1L] / config$carrying_capacity
      )
      y <- donor_lineage[index]
      k1 <- rate_1 * y
      k2 <- rate_2 * (y + dt * k1 / 2)
      k3 <- rate_2 * (y + dt * k2 / 2)
      k4 <- rate_4 * (y + dt * k3)
      donor_lineage[index + 1L] <- max(
        0,
        y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
      )
    }
    for (recipient in seq_along(backgrounds)) {
      if (recipient == donor) next
      hazard <- parameters[["h"]] * config$omega[recipient, donor] *
        s_matrix[, recipient] * donor_lineage / config$carrying_capacity
      lambda[recipient, donor] <- trapz_integral(time, hazard) / seed_size
    }
  }
  establishment <- config$within_establishment[backgrounds]
  next_generation <- 1 - exp(-sweep(lambda, 1L, establishment, `*`))
  diag(next_generation) <- 0
  reproduction <- max(Mod(eigen(next_generation, only.values = TRUE)$values))
  list(
    lambda = lambda,
    matrix = next_generation,
    reproduction = reproduction,
    seed_size = seed_size,
    horizon = horizon
  )
}

#' Calculate reference and effective within-host reproduction numbers
#'
#' @param parameters Named `h`, `gamma`, and `delta` parameters.
#' @param config Multiscale configuration.
#' @param antibiotic_end Current treatment end time.
#'
#' @return A list containing reference and current next-generation results.
#' @export
within_host_reproduction_metrics <- function(
    parameters,
    config,
    antibiotic_end = config$antibiotic_end
) {
  list(
    r0 = within_host_next_generation(
      parameters,
      config,
      antibiotic_end = config$antibiotic_start
    ),
    re = within_host_next_generation(
      parameters,
      config,
      antibiotic_end = antibiotic_end
    )
  )
}

#' Calculate one carrier's trajectory-integrated between-host reproduction
#'
#' @param trajectory Within-host trajectory.
#' @param contact_rate Daily contact rate.
#' @param susceptible_fraction Current susceptible fraction.
#' @param config Multiscale configuration.
#'
#' @return A named vector containing the integrated reproduction number,
#'   product-form comparator, carriage duration, and effective contact risk.
#' @export
between_host_reproduction <- function(
    trajectory,
    contact_rate,
    susceptible_fraction,
    config
) {
  resistant_columns <- paste0("R_", config$backgrounds$background)
  susceptible_columns <- paste0("S_", config$backgrounds$background)
  resistant <- rowSums(trajectory[resistant_columns])
  total <- resistant + rowSums(trajectory[susceptible_columns])
  fraction <- resistant / pmax(total, 1e-12)
  above <- resistant >= config$carriage_detection
  duration <- if (any(above)) max(trajectory$time[above]) else 0
  keep <- trajectory$time <= duration
  acquisition_risk <- 1 - exp(-config$inoculum_kappa * fraction)
  integrated_risk <- trapz_integral(
    trajectory$time[keep],
    acquisition_risk[keep]
  )
  average_fraction <- if (duration > 0) {
    trapz_integral(trajectory$time[keep], fraction[keep]) / duration
  } else {
    0
  }
  effective_contact_risk <- if (duration > 0) integrated_risk / duration else 0
  multiplier <- contact_rate * susceptible_fraction *
    config$between_establishment
  c(
    reproduction = multiplier * integrated_risk,
    product_comparator = multiplier * duration *
      (1 - exp(-config$inoculum_kappa * average_fraction)),
    carriage_duration = duration,
    effective_contact_risk = effective_contact_risk,
    average_resistant_fraction = average_fraction
  )
}

#' Calculate reference and effective between-host reproduction numbers
#'
#' @param parameters Named `h`, `gamma`, and `delta` parameters.
#' @param contact_rate Daily contact rate.
#' @param susceptible_fraction Current susceptible fraction.
#' @param config Multiscale configuration.
#' @param antibiotic_end Current treatment end time.
#'
#' @return A list containing reference trajectory, current trajectory, and
#'   their reproduction summaries.
#' @export
between_host_reproduction_metrics <- function(
    parameters,
    contact_rate,
    susceptible_fraction,
    config,
    antibiotic_end = config$antibiotic_end
) {
  reference <- simulate_within_host(
    parameters,
    config,
    antibiotic_end = config$antibiotic_start
  )
  current <- simulate_within_host(
    parameters,
    config,
    antibiotic_end = antibiotic_end
  )
  list(
    r0 = between_host_reproduction(reference, contact_rate, 1, config),
    re = between_host_reproduction(
      current,
      contact_rate,
      susceptible_fraction,
      config
    ),
    reference_trajectory = reference,
    current_trajectory = current
  )
}

#' Calculate all multiscale reproduction targets for one parameter set
#'
#' @param parameters Named `h`, `gamma`, and `delta` parameters.
#' @param contact_rate Daily contact rate.
#' @param susceptible_fraction Current susceptible fraction.
#' @param config Multiscale configuration.
#' @param antibiotic_end Current treatment end time.
#'
#' @return A list containing within-host details and a one-row target table.
#' @export
multiscale_reproduction_targets <- function(
    parameters,
    contact_rate,
    susceptible_fraction,
    config,
    antibiotic_end = config$antibiotic_end
) {
  within <- within_host_reproduction_metrics(parameters, config, antibiotic_end)
  between <- between_host_reproduction_metrics(
    parameters,
    contact_rate,
    susceptible_fraction,
    config,
    antibiotic_end
  )
  list(
    targets = data.frame(
      r0_within = within$r0$reproduction,
      re_within = within$re$reproduction,
      r0_between = unname(between$r0[["reproduction"]]),
      re_between = unname(between$re[["reproduction"]]),
      product_between = unname(between$re[["product_comparator"]]),
      carriage_duration = unname(between$re[["carriage_duration"]]),
      effective_contact_risk = unname(between$re[["effective_contact_risk"]])
    ),
    within = within,
    between = between
  )
}
