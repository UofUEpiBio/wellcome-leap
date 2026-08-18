# AMPLIFY Computational Prototype


This repository contains a toy computational prototype for the AMPLIFY
Wellcome Leap proposal. It connects a mechanistic agent-based
disease-transmission simulation to a missing-data-aware machine-learning
surrogate.

The prototype:

1.  simulates two antibiotic-resistant organism scenarios with
    `epiworldR::ModelSEIRCONN()`;
2.  learns individual realized secondary infections, `R_i`, from an
    agent feature `X` and organism scenario using native R `torch`; and
3.  provides one prediction wrapper that works with both inputs, `X`
    only, or scenario only.

The transmission probability is
`plogis(alpha + beta * X_i + delta * S_s)`, where `X_i` is standard
normal and `S_s` distinguishes the organism scenarios.

## Prototype workflow

``` mermaid
flowchart LR
  A["Calibrated SEIR parameters"] --> B["Simulate 5,000 populations per scenario"]
  B --> C["Agent outcomes: X, scenario, and individual Ri"]
  C --> D["70/15/15 run-level split"]
  D --> E["Random modality dropout during training"]
  E --> E1["X and scenario"]
  E --> E2["X only"]
  E --> E3["Scenario only"]
  E1 --> F["Train and validate R torch model"]
  E2 --> F
  E3 --> F
  F --> G["Save pretrained weights and preprocessing metadata"]
  G --> H["Predict with both inputs or either input alone"]
```

## Installation

Install the two runtime dependencies and the native LibTorch runtime:

``` r
install.packages(c("epiworldR", "torch"))
torch::install_torch()
```

The repository also provides:

``` sh
Rscript scripts/00_install_dependencies.R
```

## Simulation experiments

The inspectable [simulation
experiment](reports/simulation_experiment.md) is generated from Quarto
and reports the overall `R_i` distribution by scenario. The production
workflow is:

``` sh
Rscript scripts/01_smoke_test.R
Rscript scripts/02_calibrate.R 200 8
Rscript scripts/03_run_production.R 5000 8 100
```

Generated simulation data are written under ignored `data/derived/`
paths.

## ML experiment

The [masked ML experiment](reports/ml_experiment.md) documents the
70/15/15 run-level split, modality dropout, training and validation
history, and test RMSE for all three input patterns. Rendering it also
creates ignored local weights and metadata used below:

``` sh
quarto render reports/ml_experiment.qmd --to gfm
```

## Prediction examples

This example simulates scenario 1 (`lower`) with 5,000 agents whose
feature is fixed at `X = 0.1`. Fixing `X` gives enough infected
observations at exactly that value to calculate an empirical individual
reproduction count. The example uses `epiworldR::run_multiple()` for 100
replicate populations and up to eight threads. It removes the
`source = -1` pseudo-source row before calculating the empirical value.

``` r
library(torch)
source("config/simulation.R")
source("R/build_model.R")
source("R/extract_outcomes.R")
source("R/masking.R")
source("R/torch_model.R")
source("R/fit_models.R")
source("R/predict.R")

config <- readRDS("data/derived/calibration.rds")$config
x_value <- 0.1
scenario_value <- "lower"
scenario_1 <- build_seirconn_model(
  config,
  scenario = scenario_value,
  x        = rep(x_value, config$n_agents)
)
saver <- epiworldR::make_saver("reproductive")
epiworldR::run_multiple(
  scenario_1$model,
  ndays    = config$max_days,
  nsims    = 100L,
  seed     = config$base_seed,
  saver    = saver,
  verbose  = FALSE,
  nthreads = min(8L, parallel::detectCores())
)
multiple_results <- epiworldR::run_multiple_get_results(
  scenario_1$model,
  nthreads = min(8L, parallel::detectCores())
)
individual_rt <- exclude_seed_pseudo_source(multiple_results$reproductive)
empirical_mean_rt <- mean(individual_rt$rt)

model <- load_masked_model(
  weights_path  = "artifacts/masked_model.pt",
  metadata_path = "artifacts/masked_model_metadata.rds"
)
```

The same pretrained model is evaluated with `X` only, scenario only, and
both inputs. A missing input is marginalized through the corresponding
modality-dropout state learned during training. The `run_multiple()`
empirical value is the cumulative secondary-case count observed through
day 60; unlike the production extractor, its saved reproductive table
does not expose per-agent terminal states for censoring. The empirical
column is repeated to make each prediction directly comparable with the
controlled simulation.

``` r
prediction_table <- rbind(
  predict_secondary_cases(
    model,
    x = x_value
  ),
  predict_secondary_cases(
    model,
    scenario = scenario_value
  ),
  predict_secondary_cases(
    model,
    x        = x_value,
    scenario = scenario_value
  )
)
prediction_table$infected_observations <- nrow(individual_rt)
prediction_table$empirical_mean_rt <- empirical_mean_rt
prediction_table <- prediction_table[c(
  "observation_pattern",
  "x",
  "scenario",
  "infected_observations",
  "empirical_mean_rt",
  "predicted_secondary_cases"
)]
knitr::kable(
  prediction_table,
  digits = 3,
  col.names = c(
    "Model inputs",
    "X",
    "Scenario",
    "Infected observations",
    "Empirical mean Ri at day 60",
    "Predicted mean Ri"
  )
)
```

| Model inputs | X | Scenario | Infected observations | Empirical mean Ri at day 60 | Predicted mean Ri |
|:---|---:|:---|---:|---:|---:|
| x_only | 0.1 | NA | 137930 | 0.964 | 0.916 |
| scenario_only | NA | lower | 137930 | 0.964 | 1.091 |
| both | 0.1 | lower | 137930 | 0.964 | 0.917 |

## Repository documents

- [FullProposal.md](FullProposal.md): grant narrative and scientific
  context.
- [PROTOTYPE_PLAN.md](PROTOTYPE_PLAN.md): modeling, calibration, and ML
  design.
- [AGENTS.md](AGENTS.md): durable instructions for AI contributors.
