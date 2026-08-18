#' Predict an individual's expected secondary infections
#'
#' @param object Fitted or loaded masked-model bundle.
#' @param x Optional numeric vector of agent features in `[0, 1]`.
#' @param scenario Optional character vector containing `lower` or `higher`.
#'
#' @return A data frame containing supplied inputs, observation pattern, and
#'   predicted expected secondary infections.
#' @export
predict_secondary_cases <- function(
    object,
    x        = NULL,
    scenario = NULL
) {
  input <- make_inference_matrix(
    x = x,
    scenario = scenario,
    preprocessor = object$metadata$preprocessor
  )
  prediction <- predict_rate_matrix(object$model, input)
  n <- nrow(input)
  x_output <- if (is.null(x)) rep(NA_real_, n) else rep(as.numeric(x), length.out = n)
  scenario_output <- if (is.null(scenario)) {
    rep(NA_character_, n)
  } else {
    rep(as.character(scenario), length.out = n)
  }
  x_output[input[, "x_observed"] == 0] <- NA_real_
  scenario_output[input[, "scenario_observed"] == 0] <- NA_character_
  pattern <- ifelse(
    input[, "x_observed"] == 1 & input[, "scenario_observed"] == 1,
    "both",
    ifelse(input[, "x_observed"] == 1, "x_only", "scenario_only")
  )
  result <- data.frame(
    x = x_output,
    scenario = scenario_output,
    observation_pattern = pattern,
    predicted_secondary_cases = prediction,
    stringsAsFactors = FALSE
  )
  rownames(result) <- NULL
  result
}
