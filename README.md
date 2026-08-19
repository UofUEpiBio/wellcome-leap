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
flowchart TB
  subgraph calibration["Calibration, scripts/02_calibrate.R"]
    K1["Analytic intercepts for target R0 of 0.5 and 4.0"]
    K2["Bounded search over 5 candidate intercepts per scenario,<br/>each scored on the eligible mean Ri of 200 replicates"]
    K3["Calibrated alpha and delta_scenario"]
    K1 --> K2 --> K3
  end

  subgraph production["Production simulation, scripts/03_run_production.R"]
    P1["1,000 ModelSEIRCONN replicates per scenario,<br/>10,000 agents, 120 days"]
    P2["One record per infected agent:<br/>X, scenario, infection day, realized Ri"]
    P3["Analysis window: infected by day 60, at least 30 days left<br/>to transmit, at least 90% susceptible, follow-up complete"]
    P1 --> P2 --> P3
  end

  subgraph surrogate["Surrogate training, reports/ml_experiment.qmd"]
    T1["70/15/15 split over run_id, applied within each scenario"]
    T2["Training runs"]
    T3["Validation runs"]
    T4["Test runs"]
    T5["Preprocessor fit on training rows only,<br/>then applied to every partition"]
    T6["masked_count_net trained with Adam on Poisson NLL,<br/>modality dropout resampled per row in every batch"]
    T7["Validation loss averaged across scenarios, then combined<br/>across the three patterns by the dropout probabilities"]
    T8["Early stopping restores the best weights"]
    T9["Test RMSE per observation pattern"]
    T1 --> T2 --> T5 --> T6
    T1 --> T3 --> T7
    T1 --> T4 --> T9
    T6 --> T7 --> T8 --> T9
  end

  K3 --> P1
  P3 --> T1
  T8 --> S1["Save weights as .pt and preprocessing metadata as .rds"]
  S1 --> S2["predict_secondary_cases with both inputs, X only, or scenario only"]
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
Rscript scripts/03_run_production.R 1000 8 100
```

Generated simulation data are written under ignored `data/derived/`
paths.

To discard all generated simulation data, fitted model artifacts, and
Quarto caches before a fresh run:

``` sh
Rscript scripts/00_clean_generated.R
```

## ML experiment

The [masked ML experiment](reports/ml_experiment.md) documents the
70/15/15 run-level split, modality dropout, training and validation
history, and test RMSE for all three input patterns. Rendering it also
creates ignored local weights and metadata used below:

``` sh
quarto render reports/ml_experiment.qmd --to gfm
```

## Surrogate architecture

R `torch` has no built-in renderer that turns a module into a picture,
and there is no R equivalent of `torchviz` or `torchview`. What it does
generate automatically is the module summary printed below and, via
`torch::jit_trace()`, a TorchScript computation graph. The diagram is
therefore drawn by hand, and the chunk underneath it asserts that the
saved weights still have the shape the diagram claims.

``` mermaid
flowchart TB
  IN["Input row, 4 values<br/>x: standardized X, zeroed when X is missing<br/>scenario: 0 for lower or 1 for higher, zeroed when missing<br/>x_observed: 0 or 1<br/>scenario_observed: 0 or 1"]
  MD{"Modality dropout, training only<br/>both, x_only, or scenario_only at 1/3 each<br/>resampled independently for every row"}
  PASS["Row used exactly as built<br/>dropout disabled"]
  H1["hidden_1, nn_linear from 4 to 32"]
  A1["ReLU"]
  H2["hidden_2, nn_linear from 32 to 16"]
  A2["ReLU"]
  OU["output, nn_linear from 16 to 1"]
  SP["softplus, then add 1e-6"]
  LM["lambda, the expected secondary cases, strictly positive"]
  LO["Poisson negative log-likelihood with log_input = FALSE,<br/>so the output is the rate itself"]

  IN -->|"training"| MD
  IN -->|"validation and inference"| PASS
  MD --> H1
  PASS --> H1
  H1 --> A1 --> H2 --> A2 --> OU --> SP --> LM
  LM -.->|"training and validation"| LO
```

Dropout is what makes one set of weights serve all three observation
patterns. During training it hides a modality at random, so the network
must learn to read the two indicator columns and fall back on whichever
input survives. At inference nothing is hidden at random:
`predict_secondary_cases()` builds the row for the pattern the caller
actually supplied.

``` r
library(torch)
source("config/simulation.R")
source("R/build_model.R")
source("R/extract_outcomes.R")
source("R/masking.R")
source("R/torch_model.R")
source("R/fit_models.R")
source("R/predict.R")

model <- load_masked_model(
  weights_path  = "artifacts/masked_model.pt",
  metadata_path = "artifacts/masked_model_metadata.rds"
)
architecture <- model$metadata$architecture
stopifnot(
  architecture$input_dim == 4L,
  architecture$hidden_dim_1 == 32L,
  architecture$hidden_dim_2 == 16L
)
model$model
```

    An `nn_module` containing 705 parameters.

    ── Modules ─────────────────────────────────────────────────────────────────────
    • hidden_1: <nn_linear> #160 parameters
    • hidden_2: <nn_linear> #528 parameters
    • output: <nn_linear> #17 parameters

The learned tensors, read back from the saved weights:

``` r
parameters <- model$model$parameters
knitr::kable(
  data.frame(
    tensor = names(parameters),
    shape = vapply(
      parameters,
      function(tensor) paste(dim(tensor), collapse = " x "),
      character(1)
    ),
    values = vapply(parameters, function(tensor) prod(dim(tensor)), numeric(1)),
    row.names = NULL
  ),
  col.names = c("Tensor", "Shape", "Values")
)
```

| Tensor          | Shape   | Values |
|:----------------|:--------|-------:|
| hidden_1.weight | 32 x 4  |    128 |
| hidden_1.bias   | 32      |     32 |
| hidden_2.weight | 16 x 32 |    512 |
| hidden_2.bias   | 16      |     16 |
| output.weight   | 1 x 16  |     16 |
| output.bias     | 1       |      1 |

## Analysis window for individual reproduction numbers

An agent’s realized `R_i` measures that agent’s own transmissibility
only when two conditions hold. The agent must have been given a full
opportunity to transmit, and the agent must have been infected before
susceptible depletion caps how many people anyone can still infect.
Averaging offspring counts over *every* infected agent violates both
conditions and is bounded by the final-size identity: the sum of all
`R_i` values equals the number of non-seed infections, so the average
over a completed epidemic is always slightly below one no matter how
transmissible the organism is.

Every simulation therefore runs for 120 days, and only agents infected
on or before day 60 enter any analysis. Those agents always have at
least 60 remaining days to transmit, comfortably above the 30-day
minimum transmission window the configuration enforces. Agents infected
after the population drops below 90% susceptible are also excluded,
because their offspring counts are limited by the depleted susceptible
pool rather than by `X` and the organism.

## Prediction examples

This example simulates scenario 1 (`lower`) twice with 10,000 agents
whose feature is fixed first at `X = -2` and then at `X = 2`. Fixing
each population’s feature gives enough infected observations at those
exact values to calculate empirical individual reproduction counts. Each
example uses `epiworldR::run_multiple()` for 100 replicate populations
over the full 120-day horizon and up to eight threads. It removes the
`source = -1` pseudo-source row and then applies the analysis window
described above before calculating each empirical value.

``` r
config <- readRDS("data/derived/calibration.rds")$config
x_values <- c(-2, 2)
scenario_value <- "lower"
controlled_results <- lapply(seq_along(x_values), function(index) {
  x_value <- x_values[[index]]
  scenario_1 <- build_seirconn_model(
    config,
    scenario = scenario_value,
    x        = rep(x_value, config$n_agents)
  )
  saver <- epiworldR::make_saver("reproductive", "total_hist")
  epiworldR::run_multiple(
    scenario_1$model,
    ndays    = config$max_days,
    nsims    = 100L,
    seed     = config$base_seed + index - 1L,
    saver    = saver,
    verbose  = FALSE,
    nthreads = min(8L, parallel::detectCores())
  )
  multiple_results <- epiworldR::run_multiple_get_results(
    scenario_1$model,
    nthreads = min(8L, parallel::detectCores())
  )
  individual_rt <- exclude_seed_pseudo_source(multiple_results$reproductive)
  early_last_day <- vapply(
    split(multiple_results$total_hist, multiple_results$total_hist$sim_num),
    early_phase_last_day,
    numeric(1),
    config = config
  )
  eligible <- individual_ri_eligible(
    individual_rt$source_exposure_date,
    config,
    early_last_day = early_last_day[as.character(individual_rt$sim_num)]
  )
  data.frame(
    x                     = x_value,
    all_observations      = nrow(individual_rt),
    all_mean_ri           = mean(individual_rt$rt),
    infected_observations = sum(eligible),
    empirical_mean_ri     = mean(individual_rt$rt[eligible])
  )
})
controlled_summary <- do.call(rbind, controlled_results)
```

The contrast between the two empirical columns below shows why the
analysis window matters. Averaging over every infected agent returns a
value pinned just under one by the final-size identity, whereas
averaging over the analysis window recovers the individual
transmissibility implied by `X`:

``` r
knitr::kable(
  controlled_summary,
  digits = 3,
  col.names = c(
    "Simulated X",
    "All infected agents",
    "Mean Ri over all infected",
    "Agents in analysis window",
    "Mean Ri in analysis window"
  )
)
```

| Simulated X | All infected agents | Mean Ri over all infected | Agents in analysis window | Mean Ri in analysis window |
|---:|---:|---:|---:|---:|
| -2 | 518 | 0.035 | 518 | 0.035 |
| 2 | 489342 | 0.999 | 39347 | 1.832 |

The same pretrained model is evaluated with `X` only, scenario only, and
both inputs. A missing input is marginalized through the corresponding
modality-dropout state learned during training. The `run_multiple()`
reproductive table does not expose per-agent terminal states, so it
cannot be censored the way the production extractor is; with a 120-day
horizon and a day-60 infection cutoff, however, an eligible agent is
essentially certain to have finished its infectious period. The
empirical column is repeated to make each prediction directly comparable
with the controlled simulation.

``` r
prediction_table <- do.call(rbind, lapply(x_values, function(x_value) {
  rbind(
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
}))
prediction_table$requested_x <- rep(x_values, each = 3L)
prediction_table <- merge(
  prediction_table,
  controlled_summary,
  by.x = "requested_x",
  by.y = "x",
  sort = FALSE
)
prediction_table <- prediction_table[c(
  "requested_x",
  "observation_pattern",
  "x",
  "scenario",
  "infected_observations",
  "empirical_mean_ri",
  "predicted_secondary_cases"
)]
knitr::kable(
  prediction_table,
  digits = 3,
  col.names = c(
    "Simulated X",
    "Model inputs",
    "X",
    "Scenario",
    "Infected observations",
    "Empirical mean Ri in analysis window",
    "Predicted mean Ri"
  )
)
```

| Simulated X | Model inputs | X | Scenario | Infected observations | Empirical mean Ri in analysis window | Predicted mean Ri |
|---:|:---|---:|:---|---:|---:|---:|
| -2 | x_only | -2 | NA | 518 | 0.035 | 0.515 |
| -2 | scenario_only | NA | lower | 518 | 0.035 | 0.426 |
| -2 | both | -2 | lower | 518 | 0.035 | 0.114 |
| 2 | x_only | 2 | NA | 39347 | 1.832 | 8.307 |
| 2 | scenario_only | NA | lower | 39347 | 1.832 | 0.426 |
| 2 | both | 2 | lower | 39347 | 1.832 | 2.250 |

The `both` rows track the empirical values: the model separates a
strongly transmitting agent from a weakly transmitting one inside the
same organism scenario. The single-input rows behave as designed and
should not be read as errors. A scenario-only prediction averages over
the training distribution of `X`, so it returns the same value for both
requests. An `X`-only prediction averages over organisms, and because
the analysis window admits about a hundred eligible `higher` agents for
every eligible `lower` agent, that average sits close to the `higher`
scenario.

## Repository documents

- [FullProposal.md](FullProposal.md): grant narrative and scientific
  context.
- [PROTOTYPE_PLAN.md](PROTOTYPE_PLAN.md): modeling, calibration, and ML
  design.
- [AGENTS.md](AGENTS.md): durable instructions for AI contributors.
