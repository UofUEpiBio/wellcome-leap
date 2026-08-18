#' Mask-aware count regression network
#'
#' A compact R `torch` multilayer perceptron with a softplus count-rate output.
#'
#' @param input_dim Number of model inputs.
#' @param hidden_dim_1 Width of the first hidden layer.
#' @param hidden_dim_2 Width of the second hidden layer.
#'
#' @return An R `torch` module generator.
masked_count_net <- torch::nn_module(
  classname = "masked_count_net",
  initialize = function(
      input_dim    = 4L,
      hidden_dim_1 = 16L,
      hidden_dim_2 = 8L
  ) {
    self$hidden_1 <- torch::nn_linear(input_dim, hidden_dim_1)
    self$hidden_2 <- torch::nn_linear(hidden_dim_1, hidden_dim_2)
    self$output <- torch::nn_linear(hidden_dim_2, 1L)
  },
  forward = function(x) {
    x <- torch::nnf_relu(self$hidden_1(x))
    x <- torch::nnf_relu(self$hidden_2(x))
    torch::nnf_softplus(self$output(x)) + 1e-6
  }
)

#' Create a mask-aware count model
#'
#' @param input_dim Number of model inputs.
#' @param hidden_dim_1 Width of the first hidden layer.
#' @param hidden_dim_2 Width of the second hidden layer.
#'
#' @return An initialized R `torch` module.
#' @export
create_masked_count_model <- function(
    input_dim    = 4L,
    hidden_dim_1 = 16L,
    hidden_dim_2 = 8L
) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The torch package is required. Run scripts/00_install_dependencies.R.")
  }
  masked_count_net(
    input_dim = input_dim,
    hidden_dim_1 = hidden_dim_1,
    hidden_dim_2 = hidden_dim_2
  )
}

