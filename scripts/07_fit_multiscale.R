source("config/multiscale.R")
source("R/within_host.R")
source("R/omics.R")
source("R/reproduction.R")
source("R/multiscale.R")
source("R/multiscale_abm.R")

truth_path <- "data/derived/multiscale/truth_study.rds"
if (!file.exists(truth_path)) {
  stop("Run scripts/06_simulate_multiscale.R before fitting profiles.")
}
args <- commandArgs(trailingOnly = TRUE)
n_starts <- if (length(args)) as.integer(args[[1]]) else 3L
cores <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L
batch_size <- if (length(args) >= 3L) as.integer(args[[3]]) else 100L
config <- default_multiscale_config()
study <- readRDS(truth_path)
if (length(study$simulated) != nrow(study$design) ||
    length(study$truth) != nrow(study$design) ||
    any(vapply(study$simulated, inherits, logical(1), "try-error"))) {
  stop("The truth study is incomplete; rerun scripts/06_simulate_multiscale.R.")
}
print(table(study$design$site_id))
study$truth_targets <- combine_multiscale_targets(study$truth)
study$truth <- NULL
invisible(gc())
#' Fit one indexed profile from the compacted truth study
#'
#' @param index Integer profile position.
#'
#' @return A fitted multiscale profile bundle.
fit_profile <- function(index) {
  simulated <- study$simulated[[index]]
  fit <- fit_within_host_profile(
    simulated$observations,
    config,
    antibiotic_end = simulated$profile$antibiotic_end,
    n_starts = n_starts,
    seed = config$base_seed + 10000L + index
  )
  build_multiscale_profile_bundle(
    simulated,
    config,
    parameter_source = "fitted",
    fit = fit
  )
}
if (!is.finite(batch_size) || batch_size < 1L) {
  stop("batch_size must be a positive integer.")
}
if (!is.finite(cores) || cores < 1L) {
  stop("cores must be a positive integer.")
}
indices <- seq_along(study$simulated)
batches <- split(indices, ceiling(indices / batch_size))
checkpoint_dir <- file.path(
  "data/derived/multiscale/fit_checkpoints",
  paste0(
    "n", length(indices),
    "_seed", config$base_seed,
    "_starts", n_starts,
    "_batch", batch_size
  )
)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
#' Fit profiles in bounded, resumable batches
#'
#' @return A list of fitted profile bundles ordered like the truth-study design.
fit_batches <- function() {
  fitted <- vector("list", length(indices))
  cluster <- NULL
  if (cores > 1L && .Platform$OS.type != "windows") {
    cluster <- parallel::makeForkCluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
  }
  for (batch in batches) {
    checkpoint_path <- file.path(
      checkpoint_dir,
      sprintf("profiles_%05d_%05d.rds", min(batch), max(batch))
    )
    results <- if (file.exists(checkpoint_path)) {
      readRDS(checkpoint_path)
    } else if (is.null(cluster)) {
      lapply(batch, fit_profile)
    } else {
      parallel::parLapply(cluster, batch, fit_profile)
    }
    profile_ids <- vapply(
      results,
      function(result) result$profile$profile_id,
      character(1)
    )
    expected_ids <- study$design$profile_id[batch]
    if (!identical(profile_ids, expected_ids)) {
      stop("Fit checkpoint does not match the current truth-study profiles.")
    }
    if (!file.exists(checkpoint_path)) {
      saveRDS(results, checkpoint_path, compress = "gzip")
    }
    fitted[batch] <- results
    rm(results)
    invisible(gc(FALSE))
    if (!is.null(cluster)) {
      invisible(parallel::clusterCall(cluster, function() gc(FALSE)))
    }
    cat("Fitted", max(batch), "of", length(indices), "profiles.\n")
  }
  fitted
}
fitted <- fit_batches()
study$fitted <- fitted
study$simulated <- NULL
rm(fitted)
invisible(gc())

profile_rows <- lapply(study$fitted, function(bundle) {
  summary <- bundle$baseline$between$re
  data.frame(
    profile_id = bundle$profile$profile_id,
    effective_contact_risk = summary[["effective_contact_risk"]],
    carriage_duration = summary[["carriage_duration"]],
    contact_rate = bundle$profile$contact_rate,
    susceptible_fraction = bundle$profile$susceptible_fraction
  )
})
profile_data <- do.call(rbind, profile_rows)
abm <- run_multiscale_abm(
  profile_data,
  n_agents = 2000L,
  prevalence = 5 / 2000,
  days = 120L,
  seed = config$base_seed + 20000L
)
study$abm <- abm
saveRDS(study, "data/derived/multiscale/fitted_study.rds", compress = "gzip")
targets <- combine_multiscale_targets(study$fitted)
print(utils::head(targets, 12L))
cat(
  "Saved",
  length(study$fitted),
  "fitted profiles and",
  nrow(targets),
  "profile/scenario target rows.\n"
)
