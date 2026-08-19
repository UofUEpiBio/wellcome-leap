#' Predict count rates from a numeric input matrix
#'
#' @param model Trained R `torch` model.
#' @param x Numeric model-input matrix.
#' @param batch_size Prediction batch size.
#'
#' @return Numeric vector of nonnegative predicted count rates.
predict_rate_matrix <- function(
    model,
    x,
    batch_size = 8192L
) {
  model$eval()
  predictions <- numeric(nrow(x))
  torch::with_no_grad({
    starts <- seq.int(1L, nrow(x), by = batch_size)
    for (start in starts) {
      idx <- start:min(start + batch_size - 1L, nrow(x))
      input <- torch::torch_tensor(x[idx, , drop = FALSE], dtype = torch::torch_float())
      predictions[idx] <- as.numeric(model(input)$squeeze(2))
    }
  })
  predictions
}

#' Calculate mean Poisson loss in batches
#'
#' @param model R `torch` model.
#' @param x Numeric model-input matrix.
#' @param y Numeric count target.
#' @param batch_size Evaluation batch size.
#'
#' @return Scalar mean Poisson negative log-likelihood.
poisson_loss_matrix <- function(
    model,
    x,
    y,
    batch_size = 8192L
) {
  model$eval()
  total <- 0
  torch::with_no_grad({
    starts <- seq.int(1L, nrow(x), by = batch_size)
    for (start in starts) {
      idx <- start:min(start + batch_size - 1L, nrow(x))
      input <- torch::torch_tensor(x[idx, , drop = FALSE], dtype = torch::torch_float())
      target <- torch::torch_tensor(y[idx], dtype = torch::torch_float())$unsqueeze(2)
      prediction <- model(input)
      loss <- torch::nnf_poisson_nll_loss(
        prediction,
        target,
        log_input = FALSE,
        reduction = "sum"
      )
      total <- total + loss$item()
    }
  })
  total / length(y)
}

#' Scenario-balancing observation weights
#'
#' The analysis window admits very unequal numbers of agents per organism
#' scenario, so an unweighted likelihood lets the larger scenario decide what
#' the network predicts when the scenario is hidden. Balancing weights give
#' each scenario the same total weight, which makes an `X`-only prediction an
#' average over equally likely organisms instead of an average pinned to the
#' more numerous one.
#'
#' @param scenario Character vector with one scenario label per observation.
#'
#' @return A numeric vector of positive weights whose mean is one.
#' @export
scenario_balance_weights <- function(scenario) {
  if (!length(scenario)) stop("scenario must contain at least one label.")
  counts <- table(as.character(scenario))
  weights <- 1 / as.numeric(counts[as.character(scenario)])
  if (anyNA(weights)) stop("scenario must not contain missing labels.")
  weights / mean(weights)
}

#' Train the masked R torch count model
#'
#' The training likelihood and the early-stopping criterion both treat the two
#' scenarios as equally important. Scenarios contribute very unequal agent
#' counts, so a pooled objective would let the more numerous scenario decide
#' both the fitted rates and the stopping epoch. Training rows are therefore
#' weighted by [scenario_balance_weights()], and the validation Poisson loss is
#' averaged over scenarios before it is averaged over modality states.
#'
#' @param training_data Agent-level training data.
#' @param validation_data Agent-level validation data.
#' @param preprocessor Training-derived preprocessing metadata.
#' @param hidden_dim_1 Width of the first hidden layer.
#' @param hidden_dim_2 Width of the second hidden layer.
#' @param learning_rate Adam learning rate.
#' @param batch_size Training batch size.
#' @param max_epochs Maximum number of training epochs.
#' @param patience Early-stopping patience in epochs.
#' @param seed Integer R and torch seed.
#' @param modality_dropout_probabilities Named probabilities for retaining both
#'   modalities, only `X`, or only scenario during each training observation.
#' @param balance_scenarios Whether to weight training rows so that every
#'   scenario contributes the same total weight to the likelihood.
#'
#' @return A list containing the best model, training history, and architecture.
#' @export
train_masked_torch <- function(
    training_data,
    validation_data,
    preprocessor,
    hidden_dim_1                   = 16L,
    hidden_dim_2                   = 8L,
    learning_rate                  = 0.01,
    batch_size                     = 2048L,
    max_epochs                     = 100L,
    patience                       = 10L,
    seed                           = 20260818L,
    modality_dropout_probabilities = default_modality_dropout_probabilities(),
    balance_scenarios              = TRUE
) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
  training <- list(
    x = make_masked_matrix(training_data, "both", preprocessor),
    y = as.numeric(training_data$secondary_cases),
    weight = if (balance_scenarios) {
      scenario_balance_weights(training_data$scenario)
    } else {
      rep(1, nrow(training_data))
    }
  )
  model <- create_masked_count_model(
    input_dim = ncol(training$x),
    hidden_dim_1 = hidden_dim_1,
    hidden_dim_2 = hidden_dim_2,
    modality_dropout_probabilities = modality_dropout_probabilities
  )
  optimizer <- torch::optim_adam(model$parameters, lr = learning_rate)
  validation_groups <- split(
    seq_len(nrow(validation_data)),
    validation_data$scenario
  )

  best_loss <- Inf
  best_state <- NULL
  stale_epochs <- 0L
  history <- vector("list", max_epochs)

  for (epoch in seq_len(max_epochs)) {
    model$train()
    order <- sample.int(nrow(training$x))
    total_training_loss <- 0
    starts <- seq.int(1L, length(order), by = batch_size)
    for (start in starts) {
      idx <- order[start:min(start + batch_size - 1L, length(order))]
      input <- torch::torch_tensor(
        training$x[idx, , drop = FALSE],
        dtype = torch::torch_float()
      )
      target <- torch::torch_tensor(
        training$y[idx],
        dtype = torch::torch_float()
      )$unsqueeze(2)
      weight <- torch::torch_tensor(
        training$weight[idx],
        dtype = torch::torch_float()
      )$unsqueeze(2)
      optimizer$zero_grad()
      prediction <- model(input)
      loss <- (torch::nnf_poisson_nll_loss(
        prediction,
        target,
        log_input = FALSE,
        reduction = "none"
      ) * weight)$mean()
      loss$backward()
      optimizer$step()
      total_training_loss <- total_training_loss + loss$item() * length(idx)
    }

    validation_losses <- vapply(
      names(modality_dropout_probabilities),
      function(pattern) {
        input <- make_masked_matrix(validation_data, pattern, preprocessor)
        mean(vapply(
          validation_groups,
          function(rows) {
            poisson_loss_matrix(
              model,
              input[rows, , drop = FALSE],
              validation_data$secondary_cases[rows],
              batch_size = batch_size
            )
          },
          numeric(1)
        ))
      },
      numeric(1)
    )
    validation_loss <- sum(
      validation_losses * modality_dropout_probabilities
    )
    history[[epoch]] <- data.frame(
      epoch = epoch,
      training_loss = total_training_loss / nrow(training$x),
      validation_loss = validation_loss
    )

    if (validation_loss < best_loss - 1e-6) {
      best_loss <- validation_loss
      best_state <- lapply(model$state_dict(), function(tensor) tensor$clone())
      stale_epochs <- 0L
    } else {
      stale_epochs <- stale_epochs + 1L
    }
    if (stale_epochs >= patience) break
  }

  model$load_state_dict(best_state)
  history <- do.call(rbind, history[seq_len(epoch)])
  list(
    model = model,
    history = history,
    best_epoch = history$epoch[which.min(history$validation_loss)],
    architecture = list(
      input_dim = ncol(training$x),
      hidden_dim_1 = hidden_dim_1,
      hidden_dim_2 = hidden_dim_2,
      modality_dropout_probabilities = modality_dropout_probabilities
    ),
    balance_scenarios = balance_scenarios
  )
}

#' Fit the complete missing-data-aware model workflow
#'
#' Training uses analysis-eligible individual reproduction numbers only, as
#' selected by [filter_analysis_agents()].
#'
#' @param data Agent-level simulation data.
#' @param train_fraction Fraction of runs assigned to training.
#' @param validation_fraction Fraction of runs assigned to validation.
#' @param hidden_dim_1 Width of the first hidden layer.
#' @param hidden_dim_2 Width of the second hidden layer.
#' @param learning_rate Adam learning rate.
#' @param batch_size Training batch size.
#' @param max_epochs Maximum number of epochs.
#' @param patience Early-stopping patience.
#' @param seed Integer split and training seed.
#' @param modality_dropout_probabilities Named probabilities for the three
#'   supported modality states during training.
#' @param balance_scenarios Whether to weight training rows so that every
#'   scenario contributes the same total weight to the likelihood.
#'
#' @return A fitted-model bundle with model, metadata, split data, and history.
#' @export
fit_masked_model <- function(
    data,
    train_fraction                 = 0.70,
    validation_fraction            = 0.15,
    hidden_dim_1                   = 16L,
    hidden_dim_2                   = 8L,
    learning_rate                  = 0.01,
    batch_size                     = 2048L,
    max_epochs                     = 100L,
    patience                       = 10L,
    seed                           = 20260818L,
    modality_dropout_probabilities = default_modality_dropout_probabilities(),
    balance_scenarios              = TRUE
) {
  data <- filter_analysis_agents(data)
  split <- split_by_run(data, train_fraction, validation_fraction, seed)
  preprocessor <- fit_preprocessor(split$train)
  trained <- train_masked_torch(
    training_data = split$train,
    validation_data = split$validation,
    preprocessor = preprocessor,
    hidden_dim_1 = hidden_dim_1,
    hidden_dim_2 = hidden_dim_2,
    learning_rate = learning_rate,
    batch_size = batch_size,
    max_epochs = max_epochs,
    patience = patience,
    seed = seed,
    modality_dropout_probabilities = modality_dropout_probabilities,
    balance_scenarios = balance_scenarios
  )
  metadata <- list(
    preprocessor = preprocessor,
    architecture = trained$architecture,
    split_run_ids = split$run_ids,
    supported_patterns = c("both", "x_only", "scenario_only"),
    modality_dropout_probabilities = modality_dropout_probabilities,
    balance_scenarios = balance_scenarios,
    target = "secondary_cases",
    seed = seed
  )
  list(
    model = trained$model,
    metadata = metadata,
    history = trained$history,
    best_epoch = trained$best_epoch,
    split = split
  )
}

#' Evaluate masked prediction RMSE
#'
#' @param fit Fitted bundle returned by [fit_masked_model()].
#' @param data Optional evaluation data; defaults to the test split.
#'
#' @return A data frame of RMSE values by observation pattern and scenario.
#' @export
evaluate_masked_rmse <- function(
    fit,
    data = NULL
) {
  if (is.null(data)) data <- fit$split$test
  patterns <- fit$metadata$supported_patterns
  rows <- list()
  for (pattern in patterns) {
    input <- make_masked_matrix(data, pattern, fit$metadata$preprocessor)
    prediction <- predict_rate_matrix(fit$model, input)
    groups <- c("overall", sort(unique(as.character(data$scenario))))
    for (group in groups) {
      idx <- if (group == "overall") rep(TRUE, nrow(data)) else data$scenario == group
      rows[[length(rows) + 1L]] <- data.frame(
        pattern = pattern,
        scenario = group,
        n = sum(idx),
        rmse = sqrt(mean((prediction[idx] - data$secondary_cases[idx])^2)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Save a fitted masked model
#'
#' @param fit Fitted bundle returned by [fit_masked_model()].
#' @param weights_path Destination for generated torch weights.
#' @param metadata_path Destination for generated RDS metadata.
#'
#' @return The paths invisibly.
#' @export
save_masked_model <- function(
    fit,
    weights_path  = "artifacts/masked_model.pt",
    metadata_path = "artifacts/masked_model_metadata.rds"
) {
  dir.create(dirname(weights_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
  torch::torch_save(fit$model$state_dict(), weights_path)
  saveRDS(
    list(
      metadata = fit$metadata,
      history = fit$history,
      best_epoch = fit$best_epoch
    ),
    metadata_path
  )
  invisible(c(weights = weights_path, metadata = metadata_path))
}

#' Load a saved masked model
#'
#' @param weights_path Path to generated torch weights.
#' @param metadata_path Path to generated RDS metadata.
#'
#' @return A fitted-model bundle ready for prediction.
#' @export
load_masked_model <- function(
    weights_path  = "artifacts/masked_model.pt",
    metadata_path = "artifacts/masked_model_metadata.rds"
) {
  saved <- readRDS(metadata_path)
  architecture <- saved$metadata$architecture
  model <- create_masked_count_model(
    input_dim = architecture$input_dim,
    hidden_dim_1 = architecture$hidden_dim_1,
    hidden_dim_2 = architecture$hidden_dim_2,
    modality_dropout_probabilities = architecture$modality_dropout_probabilities
  )
  state <- torch::torch_load(weights_path, device = "cpu")
  model$load_state_dict(state)
  model$eval()
  list(
    model = model,
    metadata = saved$metadata,
    history = saved$history,
    best_epoch = saved$best_epoch
  )
}
