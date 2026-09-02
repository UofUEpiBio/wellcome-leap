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

Both scenarios are calibrated to be supercritical and only moderately
apart, near individual reproduction numbers of 1.5 and 3.0, and the
feature coefficient `beta` is 0.5. That is a deliberate design choice
for the missing-input problem. A prediction that sees only `X` has to
average over the organism, and a prediction that sees only the organism
has to average over `X`, so those two averages are informative only when
neither input dominates the other. A subcritical `lower` scenario would
also contribute almost no analysis-eligible agents, which would silently
turn every `X`-only prediction into a `higher`-scenario prediction.

## Prototype workflow

``` mermaid
flowchart TB
  subgraph calibration["Calibration, scripts/02_calibrate.R"]
    K1["Analytic intercepts for target R0 of 1.5 and 3.0"]
    K2["Bounded search over 5 candidate intercepts per scenario,<br/>each scored on the eligible mean Ri of 200 replicates"]
    K3["Calibrated alpha and delta_scenario"]
    K1 --> K2 --> K3
  end

  subgraph production["Production simulation, scripts/03_run_production.R"]
    P1["1,000 ModelSEIRCONN replicates per scenario,<br/>10,000 agents, 120 days"]
    P2["One record per infected agent:<br/>X, scenario, infection day, realized Ri"]
    P3["Analysis window: infected by day 60, at least 30 days left<br/>to transmit, at least 95% susceptible, follow-up complete"]
    P1 --> P2 --> P3
  end

  subgraph surrogate["Surrogate training, reports/ml_experiment.qmd"]
    T1["70/15/15 split over run_id, applied within each scenario"]
    T2["Training runs"]
    T3["Validation runs"]
    T4["Test runs"]
    T5["Preprocessor fit on training rows only,<br/>then applied to every partition"]
    T6["masked_count_net trained with Adam on Poisson NLL,<br/>modality dropout resampled per row in every batch,<br/>rows weighted so both scenarios count equally"]
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

## Multiscale pARG extension

The multiscale workflow adds a four-background within-host ODE, synthetic qPCR
and plasmid-host linkage observations, explicit within- and between-host
reference/effective reproduction numbers, a summary-coupled carriage ABM, and a
missing-modality emulator. It models carriage rather than clinical infection.
The final model is trained with every supported synthetic site represented and
tested on held-out host profiles; a separate leave-site-out diagnostic records
the current limitation under site shift.

``` sh
Rscript scripts/06_simulate_multiscale.R 4
Rscript scripts/07_fit_multiscale.R 3
Rscript scripts/08_train_multiscale_emulator.R
quarto render reports/multiscale_experiment.qmd --to gfm
```

Generated truth, observation, fitted-model, ABM, and emulator artifacts remain
under ignored `data/derived/` and `artifacts/` paths. The rendered
[multiscale experiment](reports/multiscale_experiment.md) separates mechanistic
fitting error from emulator error and records the scientific assumptions that
the toy implementation exposes.

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
`R_i` values equals the number of non-seed infections, so that average
is always slightly below one no matter how transmissible the organism
is.

Every simulation therefore runs for 120 days, and only agents infected
on or before day 60 enter any analysis. Those agents always have at
least 60 remaining days to transmit, comfortably above the 30-day
minimum transmission window the configuration enforces. Agents infected
after the population drops below 95% susceptible are also excluded,
because their offspring counts are limited by the depleted susceptible
pool rather than by `X` and the organism.

## Prediction examples

The examples below check the surrogate against held-out simulations. For
each organism scenario, one population of 10,000 agents draws its
feature from the same standard normal distribution used in production,
and `epiworldR::run_multiple()` runs 100 replicate epidemics over the
full 120-day horizon on up to eight threads with seeds that were not
used for training. The `source = -1` pseudo-source row is removed and
the analysis window described above is applied. One criterion cannot be
reproduced here: the `run_multiple()` reproductive table does not expose
per-agent terminal states, so held-out agents are not screened for
having finished their infectious period. With a 120-day horizon and a
day-60 infection cutoff, an eligible agent is essentially certain to
have recovered.

Keeping the feature random, rather than fixing it at one value for the
whole population, matters here: an agent’s own `X` determines its own
transmissibility, but a population-wide `X` also changes how fast the
epidemic depletes its susceptibles, which would move the empirical
benchmark for reasons that have nothing to do with the individual.
Empirical values are therefore taken from the agents whose feature falls
within 0.1 of each target.

``` r
config <- readRDS("data/derived/calibration.rds")$config
threads <- min(8L, parallel::detectCores())
x_targets <- c(-1, 1)
x_tolerance <- 0.1

#' Simulate held-out individual reproduction numbers for one scenario
#'
#' @param scenario Character scenario label.
#' @param index Positive integer used to offset the held-out seed.
#'
#' @return A data frame with one row per infected agent, containing the agent
#'   feature, realized secondary cases, and the analysis-window flag.
held_out_individual_ri <- function(scenario, index) {
  seed <- as.integer(config$base_seed + 500000L + index)
  set.seed(seed)
  bundle <- build_seirconn_model(
    config,
    scenario = scenario,
    x        = stats::rnorm(config$n_agents, config$x_mean, config$x_sd)
  )
  saver <- epiworldR::make_saver("reproductive", "total_hist")
  epiworldR::run_multiple(
    bundle$model,
    ndays    = config$max_days,
    nsims    = 100L,
    seed     = seed,
    saver    = saver,
    verbose  = FALSE,
    nthreads = threads
  )
  results <- epiworldR::run_multiple_get_results(bundle$model, nthreads = threads)
  individual <- exclude_seed_pseudo_source(results$reproductive)
  early_last_day <- vapply(
    split(results$total_hist, results$total_hist$sim_num),
    early_phase_last_day,
    numeric(1),
    config = config
  )
  data.frame(
    scenario = scenario,
    x = bundle$agent_data$x[
      match(individual$source, bundle$agent_data$agent_id)
    ],
    secondary_cases = individual$rt,
    eligible = individual_ri_eligible(
      individual$source_exposure_date,
      config,
      early_last_day = early_last_day[as.character(individual$sim_num)]
    ),
    stringsAsFactors = FALSE
  )
}

scenarios <- c("lower", "higher")
validation <- do.call(rbind, unname(Map(
  held_out_individual_ri,
  scenarios,
  seq_along(scenarios)
)))
rownames(validation) <- NULL
```

The contrast between the two columns below shows why the analysis window
matters. Averaging over every infected agent returns a value pinned just
under one by the final-size identity in both scenarios, whereas
averaging over the analysis window recovers the transmissibility the
calibration targeted:

``` r
window_summary <- do.call(rbind, lapply(
  split(validation, validation$scenario),
  function(data) {
    data.frame(
      scenario = data$scenario[[1]],
      all_observations = nrow(data),
      all_mean_ri = mean(data$secondary_cases),
      window_observations = sum(data$eligible),
      window_mean_ri = mean(data$secondary_cases[data$eligible])
    )
  }
))
rownames(window_summary) <- NULL
knitr::kable(
  window_summary,
  digits = 3,
  col.names = c(
    "Scenario",
    "All infected agents",
    "Mean Ri over all infected",
    "Agents in analysis window",
    "Mean Ri in analysis window"
  )
)
```

| Scenario | All infected agents | Mean Ri over all infected | Agents in analysis window | Mean Ri in analysis window |
|:---|---:|---:|---:|---:|
| higher | 961216 | 0.999 | 46854 | 2.995 |
| lower | 87753 | 0.994 | 11399 | 1.485 |

The empirical benchmark for each scenario and feature value is the mean
`R_i` of the analysis-window agents in that band:

``` r
empirical <- do.call(rbind, lapply(scenarios, function(scenario) {
  do.call(rbind, lapply(x_targets, function(x_value) {
    rows <- validation$eligible &
      validation$scenario == scenario &
      abs(validation$x - x_value) <= x_tolerance
    data.frame(
      requested_scenario = scenario,
      requested_x = x_value,
      window_observations = sum(rows),
      empirical_mean_ri = mean(validation$secondary_cases[rows]),
      stringsAsFactors = FALSE
    )
  }))
}))
```

The same pretrained model is then evaluated with `X` only, scenario
only, and both inputs. A missing input is marginalized through the
corresponding modality-dropout state learned during training.

``` r
prediction_table <- do.call(rbind, lapply(seq_len(nrow(empirical)), function(i) {
  x_value <- empirical$requested_x[[i]]
  scenario_value <- empirical$requested_scenario[[i]]
  predictions <- rbind(
    predict_secondary_cases(model, x = x_value),
    predict_secondary_cases(model, scenario = scenario_value),
    predict_secondary_cases(model, x = x_value, scenario = scenario_value)
  )
  cbind(empirical[rep(i, nrow(predictions)), ], predictions)
}))
prediction_table <- prediction_table[c(
  "requested_scenario",
  "requested_x",
  "observation_pattern",
  "x",
  "scenario",
  "window_observations",
  "empirical_mean_ri",
  "predicted_secondary_cases"
)]
knitr::kable(
  prediction_table,
  digits = 3,
  row.names = FALSE,
  col.names = c(
    "Simulated scenario",
    "Simulated X",
    "Model inputs",
    "X",
    "Scenario",
    "Agents in analysis window",
    "Empirical mean Ri",
    "Predicted mean Ri"
  )
)
```

| Simulated scenario | Simulated X | Model inputs | X | Scenario | Agents in analysis window | Empirical mean Ri | Predicted mean Ri |
|:---|---:|:---|---:|:---|---:|---:|---:|
| lower | -1 | x_only | -1 | NA | 556 | 0.853 | 1.341 |
| lower | -1 | scenario_only | NA | lower | 556 | 0.853 | 1.463 |
| lower | -1 | both | -1 | lower | 556 | 0.853 | 0.846 |
| lower | 1 | x_only | 1 | NA | 554 | 2.141 | 3.165 |
| lower | 1 | scenario_only | NA | lower | 554 | 2.141 | 1.463 |
| lower | 1 | both | 1 | lower | 554 | 2.141 | 2.079 |
| higher | -1 | x_only | -1 | NA | 2246 | 1.834 | 1.341 |
| higher | -1 | scenario_only | NA | higher | 2246 | 1.834 | 2.982 |
| higher | -1 | both | -1 | higher | 2246 | 1.834 | 1.853 |
| higher | 1 | x_only | 1 | NA | 2263 | 4.131 | 3.165 |
| higher | 1 | scenario_only | NA | higher | 2263 | 4.131 | 2.982 |
| higher | 1 | both | 1 | higher | 2263 | 4.131 | 4.275 |

The `both` rows track the empirical values: the model separates a
strongly transmitting agent from a weakly transmitting one, and one
organism from the other. The single-input rows are averages, and they
are now useful averages rather than a restatement of the dominant
scenario. An `X`-only prediction repeats itself down the two scenarios
because it cannot identify the organism, and it lands near the midpoint
of the two empirical values, because the two scenarios carry equal
weight in the training likelihood. A scenario-only prediction repeats
itself across feature values because it averages over the training
distribution of `X`, which puts it closest to the truth near `X = 0`,
slightly above it because the rate is convex in `X`, and off by roughly
the size of the feature effect at `X = -1` and `X = 1`. Every
single-input prediction stays within a factor of about 1.7 of the value
it is estimating, in either direction, and supplying both inputs closes
that gap to a few percent.

## Prototype application

`app/` is a static, dependency-free web application that carries both fitted
surrogates and evaluates them in the browser. Its default multiscale tool maps
any nonempty combination of quantitative omics, genomic linkage, and clinical
or epidemiological context to fitted within- and between-host reference and
effective reproduction numbers. A companion mechanism explorer evaluates the
simplified ODE-derived equations directly, so the learned estimate can be
compared with transparent intervention scenarios. Both produce point estimates;
the prototype does not add Bayesian inference or confidence intervals.

The legacy tab remains available. It takes the same two attributes as
[predict_secondary_cases()](R/predict.R) — the agent feature and the organism
scenario — requires at least one of them, and flags the returned individual
reproduction number as outbreak potential, borderline, or self-limiting.

The page relabels those two attributes for a clinical audience, one
indicator each and one to one. The agent feature is presented as a
colonisation density in log10 CFU/g stool, converted to a z-score against a
stated population mean and standard deviation, and the two organism
scenarios are presented as an ESBL-producing *E. coli* and a
carbapenem-resistant *K. pneumoniae*. That relabelling lives in
`app/site.json` and is presentation only — the fitted weights are untouched
— and the application tabulates the mapping on its "How it works" tab, so
any estimate it shows can be traced back to a row of the training data.

The application is published to GitHub Pages by
[.github/workflows/pages.yml](.github/workflows/pages.yml) on every push
to `main` that touches `app/`.

To work on it locally:

``` sh
python3 -m http.server --directory app 8000
```

Then open <http://localhost:8000>. A plain `file://` open will not work,
because the page fetches its JSON at runtime.

The weights for both browser models travel with the repository as
`app/model.json`, written from the generated `.pt` weights by:

``` sh
Rscript scripts/05_export_web_model.R
```

That export is the one deliberate exception to the rule that fitted
weights are not committed: without it the published page has no model. It
is a small file, it is regenerated by a committed script, and the script
verifies that a `torch`-free forward pass reproduces the `torch`
predictions before writing. Re-run it whenever the surrogate is
retrained.

Presentation text — the app name, tagline, attribution, partner
institutions, flag thresholds, organism labels, and the indicator mapping —
lives in `app/site.json` and can be edited without touching any
code.

## Repository documents

- [PROTOTYPE_PLAN.md](PROTOTYPE_PLAN.md): modeling, calibration, and ML
  design.
- [AGENTS.md](AGENTS.md): durable instructions for AI contributors.
- [app/README.md](app/README.md): prototype application layout and
  deployment.
