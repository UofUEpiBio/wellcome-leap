# AMPLIFY masked ML experiment


This experiment trains the native R `torch` surrogate and evaluates its
predictions under three observation patterns: both `X` and scenario, `X`
only, and scenario only. Mean absolute error (MAE) is the primary
predictive-performance statistic.

## Production simulation data

The model is trained on the completed outcomes from all 10,000
production simulations per organism scenario. Each simulation contains
1,000 agents; only infected agents enter the individual
reproduction-number dataset, including infected agents with zero
secondary infections.

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
  manifest_path = "../data/derived/manifest.rds",
  root_dir       = ".."
)
run_counts <- table(study$runs$scenario)
if (!identical(as.integer(run_counts[c("lower", "higher")]), c(10000L, 10000L))) {
  stop("The production manifest must contain 10,000 runs per scenario.")
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

| partition  |  runs | infected_agents | lower_agents | higher_agents |
|:-----------|------:|----------------:|-------------:|--------------:|
| train      | 13962 |        10347605 |      4409152 |       5938453 |
| validation |  2992 |         2213966 |       939898 |       1274068 |
| test       |  2993 |         2216968 |       945221 |       1271747 |

Each training observation is presented under all three masks. The four
network inputs are standardized `X`, encoded scenario, and the two
observation indicators. The network has hidden widths 16 and 8 and a
softplus output trained with Poisson negative log-likelihood.

## Fitting and validation history

``` r
fit$best_epoch
```

    [1] 1

``` r
knitr::kable(head(fit$history, 5), digits = 4)
```

| epoch | training_loss | validation_loss |
|------:|--------------:|----------------:|
|     1 |        0.9837 |          0.9825 |
|     2 |        0.9826 |          0.9825 |
|     3 |        0.9826 |          0.9825 |
|     4 |        0.9826 |          0.9826 |
|     5 |        0.9826 |          0.9825 |

``` r
knitr::kable(tail(fit$history, 5), digits = 4)
```

|     | epoch | training_loss | validation_loss |
|:----|------:|--------------:|----------------:|
| 2   |     2 |        0.9826 |          0.9825 |
| 3   |     3 |        0.9826 |          0.9825 |
| 4   |     4 |        0.9826 |          0.9826 |
| 5   |     5 |        0.9826 |          0.9825 |
| 6   |     6 |        0.9826 |          0.9825 |

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

| pattern       | scenario |       n |   mae |
|:--------------|:---------|--------:|------:|
| both          | overall  | 2216968 | 1.000 |
| both          | higher   | 1271747 | 1.009 |
| both          | lower    |  945221 | 0.988 |
| x_only        | overall  | 2216968 | 0.999 |
| x_only        | higher   | 1271747 | 1.009 |
| x_only        | lower    |  945221 | 0.986 |
| scenario_only | overall  | 2216968 | 0.979 |
| scenario_only | higher   | 1271747 | 0.987 |
| scenario_only | lower    |  945221 | 0.967 |

The overall comparison requested for the prediction interface is:

``` r
overall_mae <- mae[mae$scenario == "overall", c("pattern", "n", "mae")]
rownames(overall_mae) <- NULL
knitr::kable(overall_mae, digits = 3)
```

| pattern       |       n |   mae |
|:--------------|--------:|------:|
| both          | 2216968 | 1.000 |
| x_only        | 2216968 | 0.999 |
| scenario_only | 2216968 | 0.979 |

Differences between the three MAE values quantify the predictive cost of
masking each input. In a completed finite epidemic, the average
offspring count across all infected agents is structurally pulled toward
one, and unmodeled epidemic timing and susceptible depletion account for
substantial individual variation. The network is also fitted with
Poisson loss rather than directly minimizing MAE. A later extension
could add infection time or susceptible fraction, or define a separate
early-phase prediction target, while retaining the same masking
interface.

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
