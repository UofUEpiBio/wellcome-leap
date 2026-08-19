#' Extract a fitted masked model as plain R numeric structures
#'
#' Reads the state dictionary of a trained [masked_count_net()] module and
#' returns the layer weights, biases, and preprocessing constants as ordinary R
#' vectors. Weight matrices are returned in row-major order so that a consumer
#' outside R can reconstruct them without transposing.
#'
#' @param object Fitted or loaded masked-model bundle.
#'
#' @return A named list with `layers`, `preprocessor`, and `architecture`.
#' @export
extract_web_model <- function(object) {
  state <- lapply(object$model$state_dict(), function(tensor) {
    as.array(tensor$to(device = "cpu"))
  })
  layer_names <- c("hidden_1", "hidden_2", "output")
  missing <- setdiff(
    c(paste0(layer_names, ".weight"), paste0(layer_names, ".bias")),
    names(state)
  )
  if (length(missing)) {
    stop("State dictionary is missing: ", paste(missing, collapse = ", "))
  }

  layers <- lapply(layer_names, function(layer) {
    weight <- state[[paste0(layer, ".weight")]]
    bias <- as.numeric(state[[paste0(layer, ".bias")]])
    list(
      units = nrow(weight),
      inputs = ncol(weight),
      # Row-major: consecutive blocks of `inputs` values form one unit.
      weight = as.numeric(t(weight)),
      bias = bias
    )
  })
  names(layers) <- layer_names

  preprocessor <- object$metadata$preprocessor
  list(
    layers = layers,
    preprocessor = list(
      x_mean = as.numeric(preprocessor$x_mean),
      x_sd = as.numeric(preprocessor$x_sd),
      scenario_levels = as.list(preprocessor$scenario_levels),
      input_columns = as.character(preprocessor$input_columns)
    ),
    architecture = list(
      input_dim = as.integer(object$metadata$architecture$input_dim),
      hidden_dim_1 = as.integer(object$metadata$architecture$hidden_dim_1),
      hidden_dim_2 = as.integer(object$metadata$architecture$hidden_dim_2),
      activation = "relu",
      output_activation = "softplus",
      output_offset = 1e-6
    )
  )
}

#' Reimplement the masked-model forward pass in base R
#'
#' Mirrors the evaluation-mode forward pass of [masked_count_net()] without
#' `torch`. It exists so the browser implementation shipped in `app/` can be
#' checked against a dependency-free reference that is itself checked against
#' `torch`.
#'
#' @param web_model Structure returned by [extract_web_model()].
#' @param input Numeric four-column model-input matrix.
#'
#' @return Numeric vector of predicted expected secondary infections.
#' @export
predict_web_model <- function(
    web_model,
    input
) {
  dense <- function(values, layer, activation) {
    weight <- matrix(
      layer$weight,
      nrow = layer$units,
      ncol = layer$inputs,
      byrow = TRUE
    )
    output <- values %*% t(weight) + rep(layer$bias, each = nrow(values))
    switch(
      activation,
      relu = pmax(output, 0),
      softplus = log1p(exp(-abs(output))) + pmax(output, 0),
      stop("Unsupported activation: ", activation)
    )
  }

  values <- dense(input, web_model$layers$hidden_1, "relu")
  values <- dense(values, web_model$layers$hidden_2, "relu")
  values <- dense(values, web_model$layers$output, "softplus")
  as.numeric(values) + web_model$architecture$output_offset
}

#' Write the browser-ready model export
#'
#' Serializes the fitted weights and preprocessing constants as a small JSON
#' file consumed by the static application in `app/`. The export is a
#' deliberately small committed fixture; the generated `.pt` weights it is
#' derived from remain untracked.
#'
#' @param object Fitted or loaded masked-model bundle.
#' @param path Destination JSON path.
#' @param evaluation Optional data frame of test RMSE by observation pattern
#'   and scenario, as returned by [evaluate_masked_rmse()].
#' @param source_label Free-text provenance label stored in the export.
#'
#' @return The written path, invisibly.
#' @export
export_web_model <- function(
    object,
    path         = "app/model.json",
    evaluation   = NULL,
    source_label = "artifacts/masked_model.pt"
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite package is required to write the web model export.")
  }
  web_model <- extract_web_model(object)
  # I() keeps length-one vectors from being written as JSON scalars, which the
  # browser implementation would then index as if they were arrays.
  web_model$layers <- lapply(web_model$layers, function(layer) {
    layer$weight <- I(layer$weight)
    layer$bias <- I(layer$bias)
    layer
  })
  web_model$preprocessor$input_columns <- I(web_model$preprocessor$input_columns)
  web_model$generated <- list(
    source = source_label,
    date = format(Sys.Date()),
    target = "secondary_cases",
    supported_patterns = I(as.character(object$metadata$supported_patterns))
  )
  if (!is.null(evaluation)) {
    web_model$evaluation <- evaluation[
      c("pattern", "scenario", "n", "rmse")
    ]
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    web_model,
    path,
    digits = 12,
    auto_unbox = TRUE,
    pretty = 2
  )
  invisible(path)
}
