#' Construct the default within-host state
#'
#' @param config Multiscale configuration.
#' @param resistant Logical; whether to seed the configured donor background.
#' @param seed_size Optional resistant seed abundance.
#'
#' @return A named numeric state vector with `S_` and `R_` entries.
#' @export
initial_within_host_state <- function(
    config,
    resistant = TRUE,
    seed_size = config$initial_resistant
) {
  backgrounds <- config$backgrounds$background
  susceptible <- config$initial_total * config$backgrounds$initial_share
  resistant_state <- rep(0, length(backgrounds))
  if (resistant) {
    donor <- match(config$initial_donor, backgrounds)
    resistant_state[donor] <- seed_size
    susceptible[donor] <- susceptible[donor] - seed_size
  }
  stats::setNames(
    c(susceptible, resistant_state),
    c(paste0("S_", backgrounds), paste0("R_", backgrounds))
  )
}

#' Antibiotic activity at a given time
#'
#' @param time Numeric time.
#' @param antibiotic_start Treatment start time.
#' @param antibiotic_end Treatment end time; a value no greater than the start
#'   represents no treatment.
#'
#' @return One during treatment and zero otherwise.
#' @export
antibiotic_activity <- function(
    time,
    antibiotic_start,
    antibiotic_end
) {
  as.numeric(
    antibiotic_end > antibiotic_start &
      time >= antibiotic_start &
      time < antibiotic_end
  )
}

#' Evaluate the within-host ODE derivative
#'
#' Segregational loss is tied to gross division of plasmid-bearing bacteria,
#' while density dependence removes cells from both plasmid states. Thus loss
#' continues in a renewing population at carrying capacity but is not a generic
#' constant curing hazard.
#'
#' @param time Numeric time.
#' @param state Named within-host state.
#' @param parameters Named vector containing scalar `h`, `gamma`, and `delta`.
#' @param config Multiscale configuration.
#' @param antibiotic_end Treatment end time.
#'
#' @return Named derivative vector in the same order as `state`.
#' @export
within_host_derivative <- function(
    time,
    state,
    parameters,
    config,
    antibiotic_end = config$antibiotic_end
) {
  n_backgrounds <- nrow(config$backgrounds)
  susceptible <- state[seq_len(n_backgrounds)]
  resistant <- state[n_backgrounds + seq_len(n_backgrounds)]
  total <- sum(susceptible) + sum(resistant)
  growth_s <- config$backgrounds$growth_rate
  growth_r <- growth_s * (1 - config$backgrounds$fitness_cost)
  births_s <- growth_s * susceptible
  births_r <- growth_r * resistant
  crowding_s <- births_s * total / config$carrying_capacity
  crowding_r <- births_r * total / config$carrying_capacity
  transfer <- susceptible / config$carrying_capacity * as.numeric(
    parameters[["h"]] * config$omega %*% resistant
  )
  antibiotic <- antibiotic_activity(
    time,
    config$antibiotic_start,
    antibiotic_end
  )

  d_susceptible <- births_s - crowding_s -
    parameters[["delta"]] * antibiotic * susceptible - transfer +
    parameters[["gamma"]] * births_r
  d_resistant <- (1 - parameters[["gamma"]]) * births_r -
    crowding_r + transfer
  stats::setNames(c(d_susceptible, d_resistant), names(state))
}

#' Integrate an ODE with breakpoint-aware fourth-order Runge-Kutta
#'
#' Integration restarts at every requested time and treatment breakpoint, so a
#' Runge-Kutta step never straddles a discontinuity in antibiotic activity.
#'
#' @param derivative Function with arguments `time` and `state`.
#' @param initial_state Named numeric initial state.
#' @param output_times Numeric output times.
#' @param step Maximum integration step.
#' @param breakpoints Numeric discontinuity times.
#'
#' @return A data frame with time and state columns.
#' @export
rk4_piecewise <- function(
    derivative,
    initial_state,
    output_times,
    step,
    breakpoints = numeric()
) {
  if (length(output_times) < 1L || any(diff(output_times) < 0)) {
    stop("output_times must be a nondecreasing numeric vector.")
  }
  if (!is.finite(step) || step <= 0) stop("step must be positive.")
  boundaries <- sort(unique(c(output_times, breakpoints)))
  boundaries <- boundaries[
    boundaries >= min(output_times) & boundaries <= max(output_times)
  ]
  state <- as.numeric(initial_state)
  names(state) <- names(initial_state)
  current <- min(output_times)
  records <- list()
  if (current %in% output_times) {
    records[[length(records) + 1L]] <- c(time = current, state)
  }

  for (boundary in boundaries[boundaries > current]) {
    while (current < boundary - .Machine$double.eps^0.5) {
      dt <- min(step, boundary - current)
      step_end <- current + dt
      at_breakpoint <- any(abs(step_end - breakpoints) < 1e-10)
      left_limit <- if (at_breakpoint) step_end - 1e-10 else step_end
      k1 <- derivative(current, state)
      k2 <- derivative(current + dt / 2, state + dt * k1 / 2)
      k3 <- derivative(current + dt / 2, state + dt * k2 / 2)
      k4 <- derivative(left_limit, state + dt * k3)
      state <- state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
      if (any(!is.finite(state)) || any(state < -1e-8)) {
        stop("ODE integration produced an invalid state.")
      }
      state[state < 0] <- 0
      current <- current + dt
    }
    current <- boundary
    if (boundary %in% output_times) {
      records[[length(records) + 1L]] <- c(time = boundary, state)
    }
  }
  result <- as.data.frame(do.call(rbind, records), check.names = FALSE)
  rownames(result) <- NULL
  result
}

#' Simulate one within-host trajectory
#'
#' @param parameters Named vector containing `h`, `gamma`, and `delta`.
#' @param config Multiscale configuration.
#' @param initial_state Optional named initial state.
#' @param antibiotic_end Treatment end time.
#' @param end_time Simulation end time.
#' @param output_times Optional output grid.
#' @param step Integration step.
#'
#' @return A data frame containing time and all within-host states.
#' @export
simulate_within_host <- function(
    parameters,
    config,
    initial_state  = initial_within_host_state(config),
    antibiotic_end = config$antibiotic_end,
    end_time       = config$carriage_horizon,
    output_times   = NULL,
    step           = config$ode_step
) {
  required <- c("h", "gamma", "delta")
  if (!all(required %in% names(parameters))) {
    stop("parameters must contain h, gamma, and delta.")
  }
  if (is.null(output_times)) {
    output_times <- sort(unique(c(
      seq(config$ode_start, end_time, by = step),
      end_time,
      config$observation_days[config$observation_days <= end_time]
    )))
  }
  derivative <- function(
      time,
      state
  ) {
    within_host_derivative(
      time,
      state,
      parameters,
      config,
      antibiotic_end
    )
  }
  rk4_piecewise(
    derivative,
    initial_state,
    output_times,
    step,
    breakpoints = c(config$antibiotic_start, antibiotic_end)
  )
}

#' Trapezoidal integral
#'
#' @param time Increasing numeric time vector.
#' @param value Numeric values at `time`.
#'
#' @return Scalar trapezoidal integral.
#' @export
trapz_integral <- function(
    time,
    value
) {
  if (length(time) != length(value) || length(time) < 2L) return(0)
  sum(diff(time) * (value[-length(value)] + value[-1L]) / 2)
}
