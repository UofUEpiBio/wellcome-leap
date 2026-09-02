#' Create a complete truth-or-fitted profile bundle
#'
#' @param simulated List returned by [simulate_multiscale_profile()].
#' @param config Multiscale configuration.
#' @param parameter_source Either `"truth"` or `"fitted"`.
#' @param fit Optional result from [fit_within_host_profile()].
#'
#' @return A profile bundle containing parameters, trajectories, reproduction
#'   targets, interventions, observations, and diagnostics.
#' @export
build_multiscale_profile_bundle <- function(
    simulated,
    config,
    parameter_source = c("truth", "fitted"),
    fit              = NULL
) {
  parameter_source <- match.arg(parameter_source)
  profile <- simulated$profile
  truth_parameters <- unlist(profile[c("h", "gamma", "delta")], use.names = TRUE)
  if (parameter_source == "fitted") {
    if (is.null(fit)) stop("fit is required when parameter_source is fitted.")
    parameters <- fit$parameters
  } else {
    parameters <- truth_parameters
  }
  baseline <- multiscale_reproduction_targets(
    parameters,
    profile$contact_rate,
    profile$susceptible_fraction,
    config,
    profile$antibiotic_end
  )
  stewardship <- multiscale_reproduction_targets(
    parameters,
    profile$contact_rate,
    profile$susceptible_fraction,
    config,
    config$shorter_antibiotic_end
  )
  inhibited <- parameters
  inhibited[["h"]] <- inhibited[["h"]] * 0.5
  conjugation <- multiscale_reproduction_targets(
    inhibited,
    profile$contact_rate,
    profile$susceptible_fraction,
    config,
    profile$antibiotic_end
  )
  target_rows <- Map(
    function(result, intervention) {
      cbind(
        profile[c("profile_id", "site_id")],
        parameter_source = parameter_source,
        intervention = intervention,
        result$targets,
        stringsAsFactors = FALSE
      )
    },
    list(baseline, stewardship, conjugation),
    c("baseline", "shorter_antibiotic", "conjugation_inhibition")
  )
  list(
    profile = profile,
    parameter_source = parameter_source,
    truth_parameters = truth_parameters,
    parameters = parameters,
    observations = simulated$observations,
    truth_trajectory = simulated$truth,
    fit = fit,
    targets = do.call(rbind, target_rows),
    baseline = baseline,
    interventions = list(
      shorter_antibiotic = stewardship,
      conjugation_inhibition = conjugation
    )
  )
}

#' Simulate a four-site truth dataset
#'
#' @param config Multiscale configuration.
#' @param profiles_per_site Number of profiles per site.
#' @param fit_profiles Whether to fit each profile and build fitted bundles.
#' @param n_starts Number of fitting starts when `fit_profiles` is true.
#' @param cores Number of forked workers used for independent profiles on
#'   non-Windows systems.
#'
#' @return A list containing simulated profiles, truth bundles, and optionally
#'   fitted bundles.
#' @export
simulate_multiscale_dataset <- function(
    config,
    profiles_per_site = 8L,
    fit_profiles      = FALSE,
    n_starts          = 5L,
    cores             = 1L
) {
  cores <- as.integer(cores)
  if (!is.finite(cores) || cores < 1L) stop("cores must be a positive integer.")
  #' Apply one independent profile worker using an optional fork cluster
  #'
  #' @param values Values supplied one at a time to `worker`.
  #' @param worker Function applied to each value.
  #'
  #' @return A list in the same order as `values`.
  map_profiles <- function(
      values,
      worker
  ) {
    if (cores > 1L && .Platform$OS.type != "windows") {
      cluster <- parallel::makeForkCluster(cores)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      parallel::parLapply(cluster, values, worker)
    } else {
      lapply(values, worker)
    }
  }
  sites <- multiscale_site_table()$site_id
  design <- expand.grid(
    profile_number = seq_len(profiles_per_site),
    site_id = sites,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design$profile_id <- sprintf(
    "%s_%03d",
    design$site_id,
    design$profile_number
  )
  simulated <- map_profiles(seq_len(nrow(design)), function(index) {
    simulate_multiscale_profile(
      design$site_id[index],
      design$profile_id[index],
      config,
      config$base_seed + index
    )
  })
  truth <- map_profiles(simulated, function(profile) {
    build_multiscale_profile_bundle(profile, config)
  })
  fitted <- NULL
  if (fit_profiles) {
    fitted <- map_profiles(seq_along(simulated), function(index) {
      fit <- fit_within_host_profile(
        simulated[[index]]$observations,
        config,
        antibiotic_end = simulated[[index]]$profile$antibiotic_end,
        n_starts = n_starts,
        seed = config$base_seed + 10000L + index
      )
      build_multiscale_profile_bundle(
        simulated[[index]],
        config,
        parameter_source = "fitted",
        fit = fit
      )
    })
  }
  list(design = design, simulated = simulated, truth = truth, fitted = fitted)
}

#' Combine profile-bundle targets
#'
#' @param bundles List of profile bundles.
#'
#' @return One target data frame containing all profiles and interventions.
#' @export
combine_multiscale_targets <- function(
    bundles
) {
  result <- do.call(rbind, lapply(bundles, `[[`, "targets"))
  rownames(result) <- NULL
  result
}

#' Summarize observations into fixed-width emulator features
#'
#' @param bundle Multiscale profile bundle.
#'
#' @return A one-row feature data frame.
#' @export
summarize_profile_observations <- function(
    bundle
) {
  observations <- bundle$observations
  finite_or_na <- function(
      x,
      fun
  ) {
    x <- x[is.finite(x)]
    if (length(x)) fun(x) else NA_real_
  }
  link_fields <- grep("^link_", names(observations), value = TRUE)
  link_values <- unlist(observations[link_fields], use.names = FALSE)
  link_values <- link_values[!is.na(link_values)]
  profile <- bundle$profile
  data.frame(
    profile_id = profile$profile_id,
    site_id = profile$site_id,
    qpcr_baseline = finite_or_na(observations$qpcr_log10[observations$day == 0], mean),
    qpcr_peak = finite_or_na(observations$qpcr_log10, max),
    qpcr_day30 = finite_or_na(observations$qpcr_log10[observations$day == 30], mean),
    ecoli_day30 = finite_or_na(observations$ecoli_log10[observations$day == 30], mean),
    klebsiella_day30 = finite_or_na(
      observations$klebsiella_log10[observations$day == 30],
      mean
    ),
    linked_backgrounds = if (length(link_values)) sum(link_values > 0) else NA_real_,
    linkage_observations = if (length(link_values)) length(link_values) else NA_real_,
    antibiotic_days = profile$antibiotic_end,
    contact_rate = profile$contact_rate,
    susceptible_fraction = profile$susceptible_fraction,
    stringsAsFactors = FALSE
  )
}

#' Create an emulator-ready profile-by-intervention table
#'
#' @param bundles List of multiscale profile bundles.
#'
#' @return A data frame of features and mechanistic targets.
#' @export
make_multiscale_emulator_data <- function(
    bundles
) {
  rows <- lapply(bundles, function(bundle) {
    features <- summarize_profile_observations(bundle)
    targets <- bundle$targets
    features <- features[rep(1L, nrow(targets)), , drop = FALSE]
    features$intervention <- targets$intervention
    features$antibiotic_days <- ifelse(
      targets$intervention == "shorter_antibiotic",
      3,
      features$antibiotic_days
    )
    features$conjugation_multiplier <- ifelse(
      targets$intervention == "conjugation_inhibition",
      0.5,
      1
    )
    cbind(
      features,
      targets[c("r0_within", "re_within", "r0_between", "re_between")]
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
