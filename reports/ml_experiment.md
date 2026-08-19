# AMPLIFY masked ML experiment


This experiment trains the native R `torch` surrogate with modality
dropout and evaluates its predictions under three observation patterns:
both `X` and scenario, `X` only, and scenario only. Root mean squared
error (RMSE) is the primary predictive-performance statistic.

## Production simulation data

The model uses completed individual outcomes from all 1,000 production
simulations per organism scenario. Each simulation contains 10,000
agents observed for 60 days. Recovered infected agents enter the
individual reproduction-number dataset, including agents with zero
secondary infections; agents still exposed or infectious at the horizon
are censored.

``` r
library(torch)
torch_set_num_threads(min(8L, parallel::detectCores()))
source("../config/simulation.R")
source("../R/build_model.R")
source("../R/extract_outcomes.R")
source("../R/simulate.R")
source("../R/masking.R")
source("../R/torch_model.R")
source("../R/fit_models.R")

config <- default_simulation_config()
study <- load_simulation_batches(
  manifest_path = "../data/derived/production/manifest.rds",
  root_dir       = ".."
)
run_counts <- table(study$runs$scenario)
if (!identical(as.integer(run_counts[c("lower", "higher")]), c(1000L, 1000L))) {
  stop("The production manifest must contain 1,000 runs per scenario.")
}
agents <- study$agents[study$agents$outcome_complete, ]
rm(study)
invisible(gc())
```

## Run-level train, validation, and test split

Splitting by `run_id` prevents agents from the same simulated epidemic
from appearing in multiple partitions. Within each scenario, 70% of runs
are assigned to training, 15% to validation, and 15% to final testing.

``` r
fit <- fit_masked_model(
  agents,
  train_fraction      = 0.70,
  validation_fraction = 0.15,
  hidden_dim_1        = 16L,
  hidden_dim_2        = 8L,
  learning_rate       = 0.01,
  batch_size          = 131072L,
  max_epochs          = 30L,
  patience            = 5L,
  seed                = 20260818L
)

split_summary <- do.call(rbind, lapply(
  c("train", "validation", "test"),
  function(partition) {
    data <- fit$split[[partition]]
    data.frame(
      partition = partition,
      runs = length(unique(data$run_id)),
      infected_agents = nrow(data),
      lower_agents = sum(data$scenario == "lower"),
      higher_agents = sum(data$scenario == "higher")
    )
  }
))
rownames(split_summary) <- NULL
knitr::kable(split_summary)
```

| partition  | runs | infected_agents | lower_agents | higher_agents |
|:-----------|-----:|----------------:|-------------:|--------------:|
| train      | 1400 |         6570615 |         6087 |       6564528 |
| validation |  300 |         1404561 |         1329 |       1403232 |
| test       |  300 |         1414277 |         1341 |       1412936 |

During training, torch-native modality dropout independently assigns
each observation one of three equally likely states: both inputs
retained, scenario dropped, or `X` dropped. Both modalities are never
dropped together. The four network inputs are standardized `X`, encoded
scenario, and the two observation indicators. Dropout is disabled during
validation and inference, where the requested observation pattern is
applied explicitly. The network has hidden widths 16 and 8 and a
softplus output trained with Poisson negative log-likelihood.

## Fitting and validation history

``` r
fit$best_epoch
```

    [1] 11

``` r
knitr::kable(head(fit$history, 5), digits = 4)
```

| epoch | training_loss | validation_loss |
|------:|--------------:|----------------:|
|     1 |        0.9195 |          0.8891 |
|     2 |        0.8876 |          0.8886 |
|     3 |        0.8884 |          0.8886 |
|     4 |        0.8881 |          0.8887 |
|     5 |        0.8883 |          0.8888 |

``` r
knitr::kable(tail(fit$history, 5), digits = 4)
```

|     | epoch | training_loss | validation_loss |
|:----|------:|--------------:|----------------:|
| 12  |    12 |        0.8880 |          0.8889 |
| 13  |    13 |        0.8883 |          0.8886 |
| 14  |    14 |        0.8881 |          0.8886 |
| 15  |    15 |        0.8881 |          0.8886 |
| 16  |    16 |        0.8883 |          0.8886 |

``` r
matplot(
  fit$history$epoch,
  fit$history[c("training_loss", "validation_loss")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("#3182bd", "#e6550d"),
  xlab = "Epoch",
  ylab = "Mean Poisson negative log-likelihood"
)
legend(
  "topright",
  legend = c("Training", "Validation"),
  col = c("#3182bd", "#e6550d"),
  lty = 1,
  lwd = 2,
  bty = "n"
)
```

![Training and validation Poisson loss by
epoch.](ml_experiment_files/figure-commonmark/unnamed-chunk-4-1.png)

## Test-set RMSE by available inputs

``` r
rmse <- evaluate_masked_rmse(fit)
knitr::kable(rmse, digits = 3)
```

| pattern       | scenario |       n |  rmse |
|:--------------|:---------|--------:|------:|
| both          | overall  | 1414277 | 2.326 |
| both          | higher   | 1412936 | 2.327 |
| both          | lower    |    1341 | 1.040 |
| x_only        | overall  | 1414277 | 2.326 |
| x_only        | higher   | 1412936 | 2.327 |
| x_only        | lower    |    1341 | 1.171 |
| scenario_only | overall  | 1414277 | 2.398 |
| scenario_only | higher   | 1412936 | 2.399 |
| scenario_only | lower    |    1341 | 1.119 |

The overall comparison requested for the prediction interface is:

``` r
overall_rmse <- rmse[rmse$scenario == "overall", c("pattern", "n", "rmse")]
rownames(overall_rmse) <- NULL
knitr::kable(overall_rmse, digits = 3)
```

| pattern       |       n |  rmse |
|:--------------|--------:|------:|
| both          | 1414277 | 2.326 |
| x_only        | 1414277 | 2.326 |
| scenario_only | 1414277 | 2.398 |

Differences between the three RMSE values quantify the predictive cost
of masking each input. RMSE is aligned with conditional-mean prediction
and therefore with the expected-count interpretation of the Poisson
model. Unmodeled epidemic timing and susceptible depletion still account
for substantial individual variation. A later extension could add
infection time or susceptible fraction while retaining the same
modality-dropout interface.

The incomplete-input estimates are population averages over the missing
feature. In particular, an `X`-only prediction cannot identify the
organism and a scenario-only prediction averages over the training
distribution of `X`.

## Save local prediction artifacts

The following generated files are ignored by Git. They allow
`README.qmd` to execute its three prediction examples after this report
has been rendered.

``` r
save_masked_model(
  fit,
  weights_path  = "../artifacts/masked_model.pt",
  metadata_path = "../artifacts/masked_model_metadata.rds"
)
```
