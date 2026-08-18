# AMPLIFY masked ML experiment


This experiment trains the native R `torch` surrogate and evaluates its
predictions under three observation patterns: both `X` and scenario, `X`
only, and scenario only. Mean absolute error (MAE) is the primary
predictive-performance statistic.

## Data generation

The inspection experiment simulates 40 independent populations per
organism scenario. The same code operates on the larger production
batches without changing the model definition.

``` r
library(torch)
source("../config/simulation.R")
source("../R/build_model.R")
source("../R/extract_outcomes.R")
source("../R/simulate.R")
source("../R/masking.R")
source("../R/torch_model.R")
source("../R/fit_models.R")

config <- default_simulation_config()
study <- run_simulation_study(
  config,
  n_reps  = 40L,
  workers = 1L
)
agents <- study$agents[study$agents$outcome_complete, ]
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
  batch_size          = 4096L,
  max_epochs          = 80L,
  patience            = 8L,
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
| train      |   56 |           37217 |        15464 |         21753 |
| validation |   12 |            7235 |         3234 |          4001 |
| test       |   12 |            8169 |         3500 |          4669 |

Each training observation is presented under all three masks. The four
network inputs are standardized `X`, encoded scenario, and the two
observation indicators. The network has hidden widths 16 and 8 and a
softplus output trained with Poisson negative log-likelihood.

## Fitting and validation history

``` r
fit$best_epoch
```

    [1] 3

``` r
knitr::kable(head(fit$history, 5), digits = 4)
```

| epoch | training_loss | validation_loss |
|------:|--------------:|----------------:|
|     1 |        0.9921 |          0.9828 |
|     2 |        0.9828 |          0.9821 |
|     3 |        0.9825 |          0.9820 |
|     4 |        0.9826 |          0.9821 |
|     5 |        0.9823 |          0.9823 |

``` r
knitr::kable(tail(fit$history, 5), digits = 4)
```

|     | epoch | training_loss | validation_loss |
|:----|------:|--------------:|----------------:|
| 7   |     7 |        0.9822 |          0.9821 |
| 8   |     8 |        0.9823 |          0.9821 |
| 9   |     9 |        0.9823 |          0.9825 |
| 10  |    10 |        0.9823 |          0.9822 |
| 11  |    11 |        0.9823 |          0.9821 |

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

## Test-set MAE by available inputs

``` r
mae <- evaluate_masked_mae(fit)
knitr::kable(mae, digits = 3)
```

| pattern       | scenario |    n |   mae |
|:--------------|:---------|-----:|------:|
| both          | overall  | 8169 | 1.011 |
| both          | higher   | 4669 | 1.010 |
| both          | lower    | 3500 | 1.013 |
| x_only        | overall  | 8169 | 1.016 |
| x_only        | higher   | 4669 | 1.014 |
| x_only        | lower    | 3500 | 1.018 |
| scenario_only | overall  | 8169 | 0.993 |
| scenario_only | higher   | 4669 | 0.993 |
| scenario_only | lower    | 3500 | 0.994 |

The overall comparison requested for the prediction interface is:

``` r
overall_mae <- mae[mae$scenario == "overall", c("pattern", "n", "mae")]
rownames(overall_mae) <- NULL
knitr::kable(overall_mae, digits = 3)
```

| pattern       |    n |   mae |
|:--------------|-----:|------:|
| both          | 8169 | 1.011 |
| x_only        | 8169 | 1.016 |
| scenario_only | 8169 | 0.993 |

In this inspection run, scenario-only input has the smallest MAE,
although the three values are close. This is a useful diagnostic rather
than evidence that masking improves information. In a completed finite
epidemic, the average offspring count across all infected agents is
structurally pulled toward one, and unmodeled epidemic timing and
susceptible depletion account for substantial individual variation. The
network is also fitted with Poisson loss rather than directly minimizing
MAE. A later extension could add infection time or susceptible fraction,
or define a separate early-phase prediction target, while retaining the
same masking interface.

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
