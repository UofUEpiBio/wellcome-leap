#' Split agent data by simulation run
#'
#' @param data Agent-level simulation data with `run_id` and `scenario`.
#' @param train_fraction Fraction of runs assigned to training.
#' @param validation_fraction Fraction of runs assigned to validation.
#' @param seed Integer split seed.
#'
#' @return A list containing training, validation, and test data plus run IDs.
#' @export
split_by_run <- function(
    data,
    train_fraction      = 0.70,
    validation_fraction = 0.15,
    seed                = 20260818L
) {
  if (train_fraction <= 0 || validation_fraction <= 0 ||
      train_fraction + validation_fraction >= 1) {
    stop("Split fractions must be positive and sum to less than one.")
  }
  required <- c("run_id", "scenario")
  if (!all(required %in% names(data))) {
    stop("data must contain: ", paste(required, collapse = ", "))
  }

  runs <- unique(data[required])
  set.seed(seed)
  split_one_scenario <- function(scenario_runs) {
    ids <- sample(scenario_runs$run_id)
    n <- length(ids)
    if (n < 3L) stop("Each scenario requires at least three simulation runs.")
    n_train <- max(1L, floor(train_fraction * n))
    n_validation <- max(1L, floor(validation_fraction * n))
    if (n_train + n_validation >= n) n_train <- n - n_validation - 1L
    list(
      train = ids[seq_len(n_train)],
      validation = ids[n_train + seq_len(n_validation)],
      test = ids[(n_train + n_validation + 1L):n]
    )
  }
  per_scenario <- lapply(split(runs, runs$scenario), split_one_scenario)
  ids <- lapply(c("train", "validation", "test"), function(partition) {
    unlist(lapply(per_scenario, `[[`, partition), use.names = FALSE)
  })
  names(ids) <- c("train", "validation", "test")

  list(
    train = data[data$run_id %in% ids$train, , drop = FALSE],
    validation = data[data$run_id %in% ids$validation, , drop = FALSE],
    test = data[data$run_id %in% ids$test, , drop = FALSE],
    run_ids = ids
  )
}

#' Fit feature preprocessing metadata
#'
#' @param training_data Training-only agent data containing `x`.
#'
#' @return A named list containing training mean, standard deviation, and
#'   scenario encoding.
#' @export
fit_preprocessor <- function(training_data) {
  if (!"x" %in% names(training_data)) stop("training_data must contain x.")
  x_sd <- stats::sd(training_data$x)
  if (!is.finite(x_sd) || x_sd == 0) x_sd <- 1
  list(
    x_mean = mean(training_data$x),
    x_sd = x_sd,
    scenario_levels = c(lower = 0, higher = 1),
    input_columns = c("x", "scenario", "x_observed", "scenario_observed")
  )
}

#' Build a masked model-input matrix
#'
#' @param data Agent-level data containing `x` and `scenario`.
#' @param pattern One of `both`, `x_only`, or `scenario_only`.
#' @param preprocessor Training-derived preprocessing metadata.
#'
#' @return A numeric four-column model-input matrix.
#' @export
make_masked_matrix <- function(
    data,
    pattern,
    preprocessor
) {
  patterns <- c("both", "x_only", "scenario_only")
  if (!pattern %in% patterns) {
    stop("pattern must be one of: ", paste(patterns, collapse = ", "))
  }
  encoded <- unname(preprocessor$scenario_levels[as.character(data$scenario)])
  if (anyNA(encoded)) stop("Unknown scenario label.")
  x_standardized <- (data$x - preprocessor$x_mean) / preprocessor$x_sd
  x_observed <- as.numeric(pattern %in% c("both", "x_only"))
  scenario_observed <- as.numeric(pattern %in% c("both", "scenario_only"))

  cbind(
    x = x_standardized * x_observed,
    scenario = encoded * scenario_observed,
    x_observed = rep(x_observed, nrow(data)),
    scenario_observed = rep(scenario_observed, nrow(data))
  )
}

#' Augment training observations with every supported mask
#'
#' @param data Agent-level training data.
#' @param preprocessor Training-derived preprocessing metadata.
#' @param target Name of the count target column.
#'
#' @return A list containing an augmented input matrix, target vector, and mask
#'   labels.
#' @export
augment_mask_patterns <- function(
    data,
    preprocessor,
    target = "secondary_cases"
) {
  patterns <- c("both", "x_only", "scenario_only")
  inputs <- lapply(patterns, function(pattern) {
    make_masked_matrix(data, pattern, preprocessor)
  })
  list(
    x = do.call(rbind, inputs),
    y = rep(as.numeric(data[[target]]), times = length(patterns)),
    pattern = rep(patterns, each = nrow(data))
  )
}

#' Build row-wise inputs for public inference
#'
#' @param x Optional numeric vector of agent features.
#' @param scenario Optional character vector of scenario labels.
#' @param preprocessor Training-derived preprocessing metadata.
#'
#' @return A numeric four-column model-input matrix.
make_inference_matrix <- function(
    preprocessor,
    x        = NULL,
    scenario = NULL
) {
  if (is.null(x) && is.null(scenario)) {
    stop("Provide x, scenario, or both.")
  }
  lengths <- c(if (!is.null(x)) length(x), if (!is.null(scenario)) length(scenario))
  n <- max(lengths)
  recycle_or_stop <- function(value, name) {
    if (is.null(value)) return(rep(NA, n))
    if (length(value) == 1L) return(rep(value, n))
    if (length(value) != n) stop(name, " must have length 1 or ", n, ".")
    value
  }
  x <- as.numeric(recycle_or_stop(x, "x"))
  scenario <- as.character(recycle_or_stop(scenario, "scenario"))
  x_observed <- is.finite(x)
  scenario_observed <- !is.na(scenario) & nzchar(scenario)
  if (any(!x_observed & !scenario_observed)) {
    stop("Every prediction row must provide x, scenario, or both.")
  }
  if (any(x[x_observed] < 0 | x[x_observed] > 1)) stop("Observed x must be in [0, 1].")

  encoded <- rep(0, n)
  encoded[scenario_observed] <- unname(
    preprocessor$scenario_levels[scenario[scenario_observed]]
  )
  if (anyNA(encoded)) stop("Unknown scenario label.")
  x_standardized <- rep(0, n)
  x_standardized[x_observed] <- (
    x[x_observed] - preprocessor$x_mean
  ) / preprocessor$x_sd

  cbind(
    x = x_standardized,
    scenario = encoded,
    x_observed = as.numeric(x_observed),
    scenario_observed = as.numeric(scenario_observed)
  )
}
