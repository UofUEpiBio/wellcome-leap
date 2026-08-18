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
`plogis(alpha + beta * X_i + delta * S_s)`, where `X_i` is Uniform(0,
1), and `S_s` distinguishes the organism scenarios.

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
Rscript scripts/02_calibrate.R 200 4
Rscript scripts/03_run_production.R 10000 4 100
```

Generated simulation data are written under ignored `data/derived/`
paths.

## ML experiment

The [masked ML experiment](reports/ml_experiment.md) documents the
70/15/15 run-level split, training and validation history, and test MAE
for all three input patterns. Rendering it also creates ignored local
weights and metadata used below:

``` sh
quarto render reports/ml_experiment.qmd --to gfm
```

## Prediction examples

Load the public wrapper and locally generated model:

``` r
library(torch)
source("R/masking.R")
source("R/torch_model.R")
source("R/fit_models.R")
source("R/predict.R")

model <- load_masked_model(
  weights_path  = "artifacts/masked_model.pt",
  metadata_path = "artifacts/masked_model_metadata.rds"
)
```

### Example 1: `X` only

When scenario is unavailable, the result averages over the scenario mix
learned during training.

``` r
predict_secondary_cases(
  model,
  x = 0.80
)
```

        x scenario observation_pattern predicted_secondary_cases
    1 0.8     <NA>              x_only                  1.250876

### Example 2: scenario only

When `X` is unavailable, the result averages over the training
distribution of `X` for that scenario.

``` r
predict_secondary_cases(
  model,
  scenario = "higher"
)
```

       x scenario observation_pattern predicted_secondary_cases
    1 NA   higher       scenario_only                 0.9749071

### Example 3: `X` and scenario

Providing both supported inputs produces the most specific prediction.

``` r
predict_secondary_cases(
  model,
  x = 0.80,
  scenario = "higher"
)
```

        x scenario observation_pattern predicted_secondary_cases
    1 0.8   higher                both                  1.213123

The wrapper also accepts vectors and returns one prediction per row.

## Repository documents

- [FullProposal.md](FullProposal.md): grant narrative and scientific
  context.
- [PROTOTYPE_PLAN.md](PROTOTYPE_PLAN.md): modeling, calibration, and ML
  design.
- [AGENTS.md](AGENTS.md): durable instructions for AI contributors.
