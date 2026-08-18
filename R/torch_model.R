#' Default modality-dropout probabilities
#'
#' @return A named probability vector for both inputs, `X` only, and scenario
#'   only.
default_modality_dropout_probabilities <- function() {
  c(
    both = 1 / 3,
    x_only = 1 / 3,
    scenario_only = 1 / 3
  )
}

#' Mask-aware count regression network
#'
#' A compact R `torch` multilayer perceptron with a softplus count-rate output.
#'
#' @param input_dim Number of model inputs.
#' @param hidden_dim_1 Width of the first hidden layer.
#' @param hidden_dim_2 Width of the second hidden layer.
#' @param modality_dropout_probabilities Named probabilities for retaining both
#'   modalities, retaining only `X`, or retaining only scenario.
#'
#' @return An R `torch` module generator.
masked_count_net <- torch::nn_module(
  classname = "masked_count_net",
  initialize = function(
      input_dim                      = 4L,
      hidden_dim_1                   = 16L,
      hidden_dim_2                   = 8L,
      modality_dropout_probabilities = default_modality_dropout_probabilities()
  ) {
    self$hidden_1 <- torch::nn_linear(input_dim, hidden_dim_1)
    self$hidden_2 <- torch::nn_linear(hidden_dim_1, hidden_dim_2)
    self$output <- torch::nn_linear(hidden_dim_2, 1L)
    self$modality_dropout_probabilities <- modality_dropout_probabilities
  },
  forward = function(x) {
    if (self$training) {
      x <- apply_modality_dropout(
        x,
        self$modality_dropout_probabilities
      )
    }
    x <- torch::nnf_relu(self$hidden_1(x))
    x <- torch::nnf_relu(self$hidden_2(x))
    torch::nnf_softplus(self$output(x)) + 1e-6
  }
)

#' Apply categorical modality dropout to a torch input batch
#'
#' Each row independently retains both modalities, only `X`, or only scenario.
#' The corresponding observation indicators are updated with the values.
#'
#' @param x Four-column R `torch` tensor containing standardized `X`, encoded
#'   scenario, and two observation indicators.
#' @param probabilities Named probabilities for `both`, `x_only`, and
#'   `scenario_only`.
#'
#' @return A tensor with one of the three supported modality masks per row.
#' @export
apply_modality_dropout <- function(
    x,
    probabilities = default_modality_dropout_probabilities()
) {
  required <- c("both", "x_only", "scenario_only")
  if (!identical(names(probabilities), required)) {
    stop("probabilities must be named: both, x_only, scenario_only.")
  }
  if (any(!is.finite(probabilities)) || any(probabilities < 0) ||
      abs(sum(probabilities) - 1) > 1e-8) {
    stop("probabilities must be nonnegative and sum to one.")
  }
  probability_tensor <- torch::torch_tensor(
    unname(probabilities),
    dtype = torch::torch_float(),
    device = x$device
  )
  pattern <- torch::torch_multinomial(
    probability_tensor,
    num_samples = x$size(1),
    replacement = TRUE
  )
  x_observed <- (pattern != 3L)$type_as(x)
  scenario_observed <- (pattern != 2L)$type_as(x)
  torch::torch_stack(
    list(
      x[, 1] * x_observed,
      x[, 2] * scenario_observed,
      x_observed,
      scenario_observed
    ),
    dim = 2L
  )
}

#' Create a mask-aware count model
#'
#' @param input_dim Number of model inputs.
#' @param hidden_dim_1 Width of the first hidden layer.
#' @param hidden_dim_2 Width of the second hidden layer.
#' @param modality_dropout_probabilities Named probabilities for the three
#'   supported modality states during training.
#'
#' @return An initialized R `torch` module.
#' @export
create_masked_count_model <- function(
    input_dim                      = 4L,
    hidden_dim_1                   = 16L,
    hidden_dim_2                   = 8L,
    modality_dropout_probabilities = default_modality_dropout_probabilities()
) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The torch package is required. Run scripts/00_install_dependencies.R.")
  }
  masked_count_net(
    input_dim = input_dim,
    hidden_dim_1 = hidden_dim_1,
    hidden_dim_2 = hidden_dim_2,
    modality_dropout_probabilities = modality_dropout_probabilities
  )
}
