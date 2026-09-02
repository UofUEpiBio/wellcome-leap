#' Multiscale emulator modality blocks
#'
#' @return A named list mapping modality blocks to numeric feature names.
#' @export
multiscale_feature_blocks <- function() {
  list(
    quantitative = c(
      "qpcr_baseline", "qpcr_peak", "qpcr_day30",
      "ecoli_day30", "klebsiella_day30"
    ),
    genomic = c("linked_backgrounds", "linkage_observations"),
    clinical = c(
      "antibiotic_days", "contact_rate", "susceptible_fraction",
      "conjugation_multiplier", paste0("site_", letters[1:4])
    )
  )
}

#' Supported multiscale observation patterns
#'
#' @return A named list containing all seven nonempty modality combinations.
#' @export
multiscale_observation_patterns <- function() {
  list(
    all = c("quantitative", "genomic", "clinical"),
    no_quantitative = c("genomic", "clinical"),
    no_genomic = c("quantitative", "clinical"),
    no_clinical = c("quantitative", "genomic"),
    quantitative_only = "quantitative",
    genomic_only = "genomic",
    clinical_only = "clinical"
  )
}

#' Split multiscale data by synthetic site
#'
#' @param data Emulator-ready profile data.
#' @param test_site Site reserved for testing.
#' @param validation_site Site reserved for validation.
#'
#' @return A list containing training, validation, and test data.
#' @export
split_multiscale_by_site <- function(
    data,
    test_site       = "site_d",
    validation_site = "site_c"
) {
  if (test_site == validation_site) {
    stop("test_site and validation_site must differ.")
  }
  list(
    train = data[!data$site_id %in% c(test_site, validation_site), , drop = FALSE],
    validation = data[data$site_id == validation_site, , drop = FALSE],
    test = data[data$site_id == test_site, , drop = FALSE]
  )
}

#' Add fixed site indicators to emulator features
#'
#' @param data Emulator-ready data.
#'
#' @return `data` with four numeric site indicator columns.
add_multiscale_site_indicators <- function(
    data
) {
  for (site in paste0("site_", letters[1:4])) {
    data[[site]] <- as.numeric(data$site_id == site)
  }
  data
}

#' Fit multiscale feature preprocessing
#'
#' @param training_data Training-site emulator data.
#'
#' @return A list containing feature blocks and training-derived centers/scales.
#' @export
fit_multiscale_preprocessor <- function(
    training_data
) {
  training_data <- add_multiscale_site_indicators(training_data)
  blocks <- multiscale_feature_blocks()
  feature_names <- unlist(blocks, use.names = FALSE)
  center <- vapply(feature_names, function(field) {
    values <- training_data[[field]]
    values <- values[is.finite(values)]
    if (length(values)) mean(values) else 0
  }, numeric(1))
  scale <- vapply(feature_names, function(field) {
    values <- training_data[[field]]
    values <- values[is.finite(values)]
    value <- if (length(values) > 1L) stats::sd(values) else 1
    if (is.finite(value) && value > 0) value else 1
  }, numeric(1))
  list(blocks = blocks, feature_names = feature_names, center = center, scale = scale)
}

#' Build a block-masked multiscale input matrix
#'
#' @param data Emulator-ready data.
#' @param observed_blocks Character vector of retained modality blocks.
#' @param preprocessor Training-derived preprocessing metadata.
#'
#' @return A numeric matrix of standardized values and block indicators.
#' @export
make_multiscale_input_matrix <- function(
    data,
    observed_blocks,
    preprocessor
) {
  unknown <- setdiff(observed_blocks, names(preprocessor$blocks))
  if (length(unknown)) stop("Unknown modality blocks: ", paste(unknown, collapse = ", "))
  prepared <- add_multiscale_site_indicators(data)
  values <- matrix(
    0,
    nrow = nrow(prepared),
    ncol = length(preprocessor$feature_names),
    dimnames = list(NULL, preprocessor$feature_names)
  )
  indicators <- matrix(
    0,
    nrow = nrow(prepared),
    ncol = length(preprocessor$blocks),
    dimnames = list(NULL, paste0(names(preprocessor$blocks), "_observed"))
  )
  for (block in names(preprocessor$blocks)) {
    fields <- preprocessor$blocks[[block]]
    raw <- as.matrix(prepared[fields])
    row_observed <- rowSums(is.finite(raw)) > 0 & block %in% observed_blocks
    for (field in fields) {
      field_values <- prepared[[field]]
      field_values[!is.finite(field_values)] <- preprocessor$center[[field]]
      standardized <- (
        field_values - preprocessor$center[[field]]
      ) / preprocessor$scale[[field]]
      values[, field] <- standardized * row_observed
    }
    indicators[, paste0(block, "_observed")] <- as.numeric(row_observed)
  }
  cbind(values, indicators)
}

#' Compact multiscale regression network
#'
#' @return An R `torch` module generator predicting four log reproduction
#'   metrics.
multiscale_regression_net <- torch::nn_module(
  classname = "multiscale_regression_net",
  initialize = function(
      input_dim,
      hidden_dim_1 = 32L,
      hidden_dim_2 = 16L,
      output_dim   = 4L
  ) {
    self$hidden_1 <- torch::nn_linear(input_dim, hidden_dim_1)
    self$hidden_2 <- torch::nn_linear(hidden_dim_1, hidden_dim_2)
    self$output <- torch::nn_linear(hidden_dim_2, output_dim)
  },
  forward = function(x) {
    x <- torch::nnf_relu(self$hidden_1(x))
    x <- torch::nnf_relu(self$hidden_2(x))
    self$output(x)
  }
)

#' Create a multiscale emulator network
#'
#' @param input_dim Number of model inputs.
#' @param hidden_dim_1 First hidden width.
#' @param hidden_dim_2 Second hidden width.
#'
#' @return An initialized R `torch` model.
#' @export
create_multiscale_emulator <- function(
    input_dim,
    hidden_dim_1 = 32L,
    hidden_dim_2 = 16L
) {
  multiscale_regression_net(
    input_dim = input_dim,
    hidden_dim_1 = hidden_dim_1,
    hidden_dim_2 = hidden_dim_2
  )
}

#' Fit the leave-site-out multiscale emulator
#'
#' @param data Emulator-ready profile data.
#' @param test_site Held-out test site.
#' @param validation_site Held-out validation site.
#' @param hidden_dim_1 First hidden width.
#' @param hidden_dim_2 Second hidden width.
#' @param learning_rate Adam learning rate.
#' @param max_epochs Maximum epochs.
#' @param patience Early-stopping patience.
#' @param seed Random seed.
#'
#' @return A fitted model bundle with preprocessing, splits, and history.
#' @export
fit_multiscale_emulator <- function(
    data,
    test_site       = "site_d",
    validation_site = "site_c",
    hidden_dim_1    = 32L,
    hidden_dim_2    = 16L,
    learning_rate   = 0.01,
    max_epochs      = 200L,
    patience        = 20L,
    seed            = 20260901L
) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
  split <- split_multiscale_by_site(data, test_site, validation_site)
  preprocessor <- fit_multiscale_preprocessor(split$train)
  patterns <- multiscale_observation_patterns()
  target_names <- c("r0_within", "re_within", "r0_between", "re_between")
  make_augmented <- function(
      partition
  ) {
    inputs <- lapply(patterns, function(blocks) {
      make_multiscale_input_matrix(partition, blocks, preprocessor)
    })
    list(
      x = do.call(rbind, inputs),
      y = log(pmax(
        as.matrix(partition[target_names])[rep(seq_len(nrow(partition)), length(patterns)), ],
        1e-8
      ))
    )
  }
  training <- make_augmented(split$train)
  validation <- make_augmented(split$validation)
  model <- create_multiscale_emulator(
    ncol(training$x),
    hidden_dim_1,
    hidden_dim_2
  )
  optimizer <- torch::optim_adam(model$parameters, lr = learning_rate)
  train_x <- torch::torch_tensor(training$x, dtype = torch::torch_float())
  train_y <- torch::torch_tensor(training$y, dtype = torch::torch_float())
  validation_x <- torch::torch_tensor(validation$x, dtype = torch::torch_float())
  validation_y <- torch::torch_tensor(validation$y, dtype = torch::torch_float())
  best_loss <- Inf
  best_state <- NULL
  stale <- 0L
  history <- vector("list", max_epochs)
  for (epoch in seq_len(max_epochs)) {
    model$train()
    optimizer$zero_grad()
    training_loss <- torch::nnf_mse_loss(model(train_x), train_y)
    training_loss$backward()
    optimizer$step()
    model$eval()
    torch::with_no_grad({
      validation_loss <- torch::nnf_mse_loss(model(validation_x), validation_y)
    })
    value <- validation_loss$item()
    history[[epoch]] <- data.frame(
      epoch = epoch,
      training_loss = training_loss$item(),
      validation_loss = value
    )
    if (value < best_loss - 1e-6) {
      best_loss <- value
      best_state <- lapply(model$state_dict(), function(tensor) tensor$clone())
      stale <- 0L
    } else {
      stale <- stale + 1L
    }
    if (stale >= patience) break
  }
  model$load_state_dict(best_state)
  list(
    model = model,
    preprocessor = preprocessor,
    split = split,
    patterns = patterns,
    target_names = target_names,
    history = do.call(rbind, history[seq_len(epoch)]),
    architecture = list(
      input_dim = ncol(training$x),
      hidden_dim_1 = hidden_dim_1,
      hidden_dim_2 = hidden_dim_2
    )
  )
}

#' Predict multiscale reproduction metrics
#'
#' @param object Fitted multiscale emulator.
#' @param data Emulator-ready input data.
#' @param pattern Named observation pattern.
#'
#' @return A data frame of positive reproduction predictions.
#' @export
predict_multiscale_emulator <- function(
    object,
    data,
    pattern = "all"
) {
  blocks <- object$patterns[[pattern]]
  if (is.null(blocks)) stop("Unknown observation pattern: ", pattern)
  input <- make_multiscale_input_matrix(data, blocks, object$preprocessor)
  object$model$eval()
  torch::with_no_grad({
    prediction <- object$model(torch::torch_tensor(
      input,
      dtype = torch::torch_float()
    ))
  })
  result <- as.data.frame(exp(as.matrix(prediction)))
  names(result) <- object$target_names
  result
}

#' Evaluate leave-site-out emulator performance
#'
#' @param object Fitted multiscale emulator.
#' @param data Optional evaluation data; defaults to the test site.
#'
#' @return A data frame of RMSE and threshold balanced accuracy by target and
#'   observation pattern.
#' @export
evaluate_multiscale_emulator <- function(
    object,
    data = object$split$test
) {
  rows <- list()
  for (pattern in names(object$patterns)) {
    predicted <- predict_multiscale_emulator(object, data, pattern)
    for (target in object$target_names) {
      truth <- data[[target]]
      prediction <- predicted[[target]]
      class <- truth > 1
      predicted_class <- prediction > 1
      sensitivity <- if (any(class)) mean(predicted_class[class]) else NA_real_
      specificity <- if (any(!class)) mean(!predicted_class[!class]) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        pattern = pattern,
        target = target,
        n = length(truth),
        rmse = sqrt(mean((prediction - truth)^2)),
        balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE)
      )
    }
  }
  do.call(rbind, rows)
}

#' Save a fitted multiscale emulator
#'
#' @param object Fitted multiscale emulator.
#' @param weights_path Destination for generated torch weights.
#' @param metadata_path Destination for generated metadata.
#'
#' @return Paths invisibly.
#' @export
save_multiscale_emulator <- function(
    object,
    weights_path  = "artifacts/multiscale_emulator.pt",
    metadata_path = "artifacts/multiscale_emulator_metadata.rds"
) {
  dir.create(dirname(weights_path), recursive = TRUE, showWarnings = FALSE)
  torch::torch_save(object$model$state_dict(), weights_path)
  saveRDS(
    object[c(
      "preprocessor", "patterns", "target_names", "history", "architecture"
    )],
    metadata_path
  )
  invisible(c(weights = weights_path, metadata = metadata_path))
}

#' Load a fitted multiscale emulator
#'
#' @param weights_path Path to generated torch weights.
#' @param metadata_path Path to generated metadata.
#'
#' @return A model bundle ready for [predict_multiscale_emulator()].
#' @export
load_multiscale_emulator <- function(
    weights_path  = "artifacts/multiscale_emulator.pt",
    metadata_path = "artifacts/multiscale_emulator_metadata.rds"
) {
  saved <- readRDS(metadata_path)
  model <- create_multiscale_emulator(
    saved$architecture$input_dim,
    saved$architecture$hidden_dim_1,
    saved$architecture$hidden_dim_2
  )
  model$load_state_dict(torch::torch_load(weights_path, device = "cpu"))
  model$eval()
  c(list(model = model), saved)
}
