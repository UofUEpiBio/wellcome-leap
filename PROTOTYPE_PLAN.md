# AMPLIFY Toy Prototype: Simulation and Missing-Data ML Plan

## 1. Purpose

Build a small, reproducible demonstration of the computational workflow proposed in AMPLIFY:

1. use a mechanistic agent-based SEIR model to generate heterogeneous transmission outcomes;
2. calibrate scenario-level epidemic potential;
3. recover individual transmission outcomes and observable features; and
4. train a machine-learning surrogate that predicts individual transmission under incomplete feature availability.

This prototype represents the between-host and population-level portions of Thrust 1A and provides a minimal demonstration of the proposal's missing-modality ML strategy. It is not intended to represent plasmid biology or the full multiscale system.

## 2. Scope of the first prototype

### Simulation engine

- R with `epiworldR` (installed version: 0.15.1.0).
- `ModelSEIRCONN`, a fully connected stochastic SEIR model.
- Population size: 1,000 agents per simulation.
- Two transmission scenarios, each representing a different antibiotic-resistant organism.
- A calibration stage followed by 10,000 retained simulations per scenario (20,000 total).
- Reproducible, independent random-number streams and parallel execution.

### Agent-level heterogeneity

For agent \(i\), draw

\[
X_i \sim \operatorname{Uniform}(0,1).
\]

The first prototype will use the agent-feature logit function provided natively by `epiworldR`:

\[
p_{is} = \operatorname{logit}^{-1}(\alpha + \beta X_i + \delta S_s),
\]

where:

- \(p_{is}\) is agent \(i\)'s per-contact probability of transmitting in scenario \(s\);
- \(S_s=0\) for the lower-transmission organism and \(S_s=1\) for the higher-transmission organism;
- \(\alpha\) is the baseline log-odds for the lower-transmission organism;
- \(\beta>0\) controls the strength of continuous agent heterogeneity, provisionally \(\beta=1\); and
- \(\delta>0\) shifts the probability upward for the higher-transmission organism.

This retains the requested monotonic, nonlinear relationship between \(X\) and transmissibility while guaranteeing a valid probability. On the odds scale, a one-unit increase in \(X\) multiplies the transmission odds by \(e^\beta\), and changing from scenario 1 to scenario 2 multiplies the odds by \(e^\delta\).

### Feasibility with the R interface

No custom C++ model is required for this functional form. The implementation path available in `epiworldR` 0.15.1.0 is:

1. Construct the standard `ModelSEIRCONN` model.
2. Attach an agent-data matrix with columns `intercept = 1`, `x`, and `scenario_high` using `set_agents_data()`.
3. Create a callback with `virus_fun_logit(vars = c(0L, 1L, 2L), coefs = c(alpha, beta, delta), model = model)`; variable indices are zero-based.
4. Attach it to the model's virus using `set_prob_infecting_fun(get_virus(model, 0), model, vfun)`.

The relevant source path has been verified: `ModelSEIRCONN` obtains the virus from the infectious source agent, and the virus callback evaluates its probability using the agent that owns that virus. Thus, the callback reads \(X_i\) from the transmitter, not from the susceptible recipient.

A temporary in-memory feasibility experiment using 100 replicates per scenario confirmed that this route produces both a positive \(X\) effect and a scenario shift in individual secondary infections. The lower and higher pilot scenarios produced mean early \(R_i\) values of 1.39 and 1.82, respectively. These values are evidence that the mechanism works, not final calibrated estimates.

Custom C++ `epiworld` should be reserved as a fallback if the project later requires an arbitrary non-logistic function, an interaction involving both source and recipient features, time-varying agent features, or other dynamics not exposed by the R callback factories.

### Epidemiological constants

The first implementation will keep incubation, recovery, initial prevalence, contact rate, and run horizon constant across scenarios except where a calibration decision below requires otherwise. These values will be stored in one configuration file rather than embedded in analysis code.

Provisional values for an initial smoke test, subject to approval:

| Parameter | Provisional value | Rationale |
|---|---:|---|
| Population | 1,000 | Requested |
| Initial infected/exposed seeds | 1% (10 agents) | Reduces early stochastic extinction during pipeline testing |
| Mean incubation period | 3 days | Simple toy-model value |
| Daily recovery probability | 0.20 | Mean infectious duration of about 5 days |
| Simulation horizon | At least 100 days and then until extinction, subject to a safety cap | Prevents incomplete follow-up of late infections |
| \(X\) distribution | Uniform(0,1) | Requested |
| \(\beta\) | 1 | Gives an \(e\)-fold transmission-odds ratio across the range of \(X\) |

These are computational defaults, not estimates for Enterobacterales or pARG transmission.

## 3. Terminology and estimands

`epiworldR::get_reproductive_number()` returns one record per infected source, including sources with zero onward transmissions. Its `rt` column is the individual realized secondary-case count used here as \(R_i\). The transmission table produced by `get_transmissions()` identifies the corresponding source and target of every infection and will be used as an extraction cross-check. For an infected agent \(i\), the count

\[
R_i = \sum_j I(i \text{ infected } j)
\]

is the agent's **realized secondary-case count** (also called an individual offspring count or individual reproduction number).

This is not the conventional population-level time-varying reproduction number \(R_t\), which is indexed by calendar time and normally estimated from incidence and a generation-interval distribution. To avoid an ambiguous scientific claim, the prototype should call the primary target `secondary_cases` or \(R_i\). If desired, a separate population-level \(R_t\) time series can be estimated as a secondary output.

The population-level calibration target also needs a precise definition. The recommended prototype definition is the mean number of secondary infections caused by early infectious agents while susceptible depletion is negligible. This estimates an empirical early-epidemic reproduction number and is more appropriate for calibration than averaging offspring counts over an entire finite epidemic.

## 4. Experimental design

### Recommended interpretation

Run the two organism scenarios separately, with one transmissible process in each run:

- Scenario 1: lower-transmission antibiotic-resistant organism, target early-epidemic \(R_0 \approx 1.5\).
- Scenario 2: higher-transmission antibiotic-resistant organism, target early-epidemic \(R_0 \approx 2.0\).

The organism/scenario label and agent-level \(X\) will determine the transmission probability. Separate runs avoid competition between two organisms for the same susceptible population and produce a clean supervised-learning dataset.

### Redrawing agents

The recommended design redraws the 1,000 \(X_i\) values for every replicate. This makes the replicates independent population realizations and gives the ML model broad coverage of the continuous feature. A secondary sensitivity analysis can hold one population fixed across replicates to isolate epidemic stochasticity from population-composition variability.

### Initial infections

Seed agents must be selected independently of \(X\) unless preferential seeding is a scientific feature of interest. Seed infections appear in `get_transmissions()` with source `-1` and will not count as anyone's secondary infection.

### Retained agent records

The dataset will include **every infected agent**, including agents who caused zero secondary infections. Restricting the data to agents who infected at least one person would create a zero-truncated outcome, bias the relationship between \(X\) and transmission, and prevent the model from learning non-transmission.

Minimum long-form schema:

| Variable | Meaning |
|---|---|
| `run_id` | Unique replicate identifier |
| `seed` | Reproducibility seed |
| `scenario` | Lower- or higher-transmission scenario |
| `organism` | Antibiotic-resistant organism label associated with the scenario |
| `agent_id` | Agent identifier within replicate |
| `x` | Agent's continuous feature |
| `p_transmit` | Mechanistic per-contact transmission probability |
| `infection_day` | Day infected/exposed |
| `source_id` | Infector, where applicable |
| `secondary_cases` | Number of targets infected by this agent |
| `outcome_complete` | Whether the agent had complete transmission follow-up |
| `early_phase` | Whether the record belongs to the calibration window |
| `final_epidemic_size` | Replicate-level epidemic size |
| `extinct_early` | Indicator of stochastic early extinction |

We will also retain a run-level summary table and, optionally, the raw transmission-edge table. Raw state histories need not be retained for all 20,000 runs unless they are required for a later analysis.

## 5. Calibration strategy

Calibration will occur before the 10,000-run production stage.

1. Fix incubation, recovery, contact rate, initial seeding, \(\beta\), and the definition of the early epidemic window.
2. Obtain analytic starting values for \(\alpha\) and \(\delta\). Early in a fully susceptible epidemic, the approximation

   \[
   R_{0,s} \approx \frac{C}{\gamma}E[p_{is}]
   \]

   can be combined with the exact mean of the logistic probability for \(X\sim U(0,1)\):

   \[
   E\left[\operatorname{logit}^{-1}(a_s+\beta X)\right]
   = \frac{\log(1+e^{a_s+\beta})-\log(1+e^{a_s})}{\beta},
   \]

   where \(a_1=\alpha\), \(a_2=\alpha+\delta\), \(C\) is the daily contact rate, and \(\gamma\) is the daily recovery probability. With \(C=2\), \(\gamma=0.2\), and \(\beta=1\), solving this approximation gives initial values \(\alpha\approx-2.264\) and \(\delta\approx0.352\) for targets 1.5 and 2.0. These are starting values only because the simulator is discrete, stochastic, and finite.
3. For candidate values around that starting point, run a moderate batch of simulations using common random numbers across candidates to reduce Monte Carlo noise.
4. Estimate early-epidemic \(R_0\) from `get_reproductive_number()` among agents infected before a specified susceptible-depletion threshold (recommended: while at least 90% of the population remains susceptible).
5. Use stochastic root finding or a bounded grid/refinement search to approach each target. Exact equality is not required; the aim is clear separation around \(R_0\approx1.5\) and \(R_0\approx2.0\), with Monte Carlo intervals reported.
6. Confirm that \(\beta>0\) and \(\delta>0\), so transmission increases with \(X\) and is higher in scenario 2 while all other epidemiological parameters remain fixed.
7. Validate the chosen parameters in a fresh simulation batch not used during calibration.
8. Freeze the calibrated configuration before production runs.

## 6. Simulation workflow

### Phase A: executable smoke test

- Construct one scenario with 1,000 agents.
- Verify the distribution of \(X\) and calculated probabilities.
- Run a small number of replicates.
- Confirm transmission edges link valid sources and targets.
- Confirm individual secondary-case counts reproduce the number of non-seed transmission edges exactly.
- Confirm that zero-secondary-case outcomes have complete follow-up and are not artifacts of the simulation endpoint.
- Confirm identical seeds reproduce identical results.

### Phase B: calibration pilot

- Run approximately 200-1,000 replicates per candidate setting, increasing the number near the solution.
- Quantify calibration uncertainty and early-extinction frequency.
- Produce diagnostic plots of empirical \(R_0\), epidemic size, and secondary cases versus \(X\).

### Phase C: production simulation

- Run 10,000 independent replicates per scenario.
- Continue each epidemic until no exposed or infectious agents remain. If a prespecified safety cap is reached, flag incomplete individual outcomes as right-censored rather than treating them as genuine zeros.
- Parallelize at the replicate or batch level with deterministic random-number streams.
- Write results in batches to avoid holding all raw histories in memory.
- Save compact columnar data (preferably Parquet) plus CSV summaries for accessibility.
- Record R/package versions, configuration, timestamps, and seed ranges.

### Phase D: simulation validation

- Check achieved \(R_0\) and its Monte Carlo interval by scenario.
- Verify that mean transmission increases monotonically with \(X\).
- Verify scenario separation.
- Examine early extinction and final epidemic-size distributions.
- Check for probability saturation or invalid values.
- Confirm that results are stable across parallel worker counts.

### Phase E: Quarto simulation report

- Run a reproducible inspection experiment from `reports/simulation_experiment.qmd` rather than an ad hoc analysis script.
- Render with `format: gfm` and commit the resulting Markdown and figure assets, but not the underlying simulated agent data.
- Report the number of infected agents, mean, standard deviation, median, interquartile range, selected quantiles, and zero-secondary-infection fraction by scenario.
- Plot the overall distribution of individual \(R_i\) values by scenario and summarize achieved early-epidemic reproduction numbers.

## 7. Missing-data ML surrogate

### Prediction target and split

- Primary target: `secondary_cases`, a non-negative and typically overdispersed count.
- Split data by `run_id`, not by agent row, so agents from one simulated epidemic cannot appear in both training and evaluation data.
- Preserve scenario balance across train, validation, and test partitions.

### Models

Use a transparent statistical model as a scientific baseline, followed by one small neural surrogate implemented with the native R `torch` package (the R interface to PyTorch/LibTorch):

1. **Baseline:** Poisson and negative-binomial regression with \(X\), scenario, and their interaction. Pattern-specific reduced baselines (\(X\)-only and scenario-only) will show the information available under each mask.
2. **Primary ML model:** one compact masked multilayer perceptron built with R `torch`. It will accept four inputs—standardized `x`, encoded `scenario`, `x_observed`, and `scenario_observed`—and return a nonnegative expected secondary-case count through a softplus output. The initial architecture will be deliberately small (4 inputs, two hidden layers no wider than 16 units, and 1 output) and will be trained with Poisson negative log-likelihood.

The model will run on CPU by default. GPU support is not required for this two-feature prototype.

### Masking design

Generate masks during model training, not only during evaluation. Provisionally sample among the supported observation patterns:

- both \(X\) and scenario observed;
- only \(X\) observed; and
- only scenario observed.

The input will include `x`, encoded `scenario`, `x_observed`, and `scenario_observed`. After estimating preprocessing statistics from the training split only, observed \(X\) values will be standardized. Masked values will then be zero-filled; the observation indicators distinguish a masked value from a genuine standardized zero or scenario 0. The test set will be evaluated separately under all three fixed patterns. The mask probabilities should reflect intended real-world availability when those expectations are known; otherwise use equal representation for the initial demonstration.

The prototype need not support both features missing unless a grand-mean fallback is desired.

Missing inputs cannot be reconstructed from nothing. When scenario is absent, an \(X\)-only prediction will average over the scenario mix learned from the training data; it cannot distinguish the two organisms. When \(X\) is absent, a scenario-only prediction will average over the population distribution of \(X\) for that organism. Predictions should therefore become more specific and generally more accurate when both inputs are observed. This graceful degradation is the behavior the masking experiment is intended to demonstrate.

### Explicit R `torch` workflow

1. Install the R package with `install.packages("torch")`, then install its LibTorch/Lantern runtime with `torch::install_torch()`. Restart R if requested and verify with `torch::torch_is_installed()`.
2. Set the R and torch seeds before splitting and training. Use CPU tensors by default for reproducibility and portability.
3. Split by `run_id`, calculate the \(X\) mean and standard deviation from the training partition only, and store these preprocessing values with the model metadata.
4. During each training epoch, augment observations with the three supported mask patterns. Each original training row should be seen under all three patterns or under a balanced random sample of them.
5. Convert the four-column input matrix and count target to `torch_float()` tensors. Use mini-batches only if the retained training data are too large for memory.
6. Define the network with `nn_module()`, `nn_linear()`, a simple activation such as ReLU, and a final `nn_softplus()` layer. Avoid dropout and architectural complexity unless validation shows a clear need.
7. Train with `optim_adam()`, Poisson negative log-likelihood, early stopping on validation loss, and restoration of the best weights.
8. Save the network `state_dict` with `torch_save()` and save a small RDS metadata bundle containing preprocessing statistics, scenario encoding, architecture parameters, and supported mask patterns. Do not commit either generated artifact.
9. Reconstruct the network and load weights on CPU with `torch_load(..., device = "cpu")`. The public wrapper must accept `x`, `scenario`, or both; reject calls where both are missing.
10. Test serialization round trips and prediction equality before and after loading.

### Evaluation

Report performance separately for each observation pattern and scenario:

- mean absolute error (MAE), the primary predictive-performance statistic;
- root mean squared error;
- Poisson deviance or negative-binomial log score;
- calibration of predicted versus observed mean counts;
- performance for zero versus nonzero secondary transmission; and
- uncertainty intervals or bootstrap intervals for aggregate metrics.

Use a run-level 70% training, 15% validation, and 15% test split. The ML Quarto report must show the split sizes, training/validation loss history, stopping epoch, and test MAE for: both \(X\) and scenario, \(X\) only, and scenario only. Also compare predictions with the known conditional simulation relationship. Because the simulator is the ground truth in this toy study, recovery of the correct monotonic effect and sensible degradation under masking are more important than maximizing a single accuracy score.

## 8. Proposed repository structure

```text
wellcome-leap/
├── AGENTS.md
├── FullProposal.md
├── PROTOTYPE_PLAN.md
├── README.qmd
├── README.md
├── config/
│   └── simulation.R
├── R/
│   ├── build_model.R
│   ├── calibrate.R
│   ├── simulate.R
│   ├── extract_outcomes.R
│   ├── masking.R
│   ├── torch_model.R
│   ├── fit_models.R
│   └── predict.R
├── scripts/
│   ├── 00_install_dependencies.R
│   ├── 01_smoke_test.R
│   ├── 02_calibrate.R
│   ├── 03_run_production.R
│   └── 04_fit_and_evaluate.R
├── tests/testthat/
├── data/
│   └── README.md
├── artifacts/
│   └── README.md
└── reports/
    ├── simulation_experiment.qmd
    ├── ml_experiment.qmd
    └── prediction_example.qmd
```

Generated simulation data should be excluded from version control unless a small example dataset is deliberately committed.

## 9. Deliverables

1. Reproducible R environment and configuration.
2. Tested `ModelSEIRCONN` builder with continuous agent-specific transmissibility.
3. Calibration script and calibration report.
4. Parallel simulation pipeline for 20,000 retained runs.
5. Agent-level, run-level, and optional edge-level datasets with data dictionaries.
6. Missingness-augmented statistical and ML models.
7. Evaluation report with simulation diagnostics, masking-stratified predictive performance, and limitations.
8. A small example run that can execute quickly on a laptop.

## 10. Implementation gates

The scientific design decisions in Gate 1 have been resolved. Implementation can begin after confirmation of the remaining operational question about the production compute environment.

### Gate 1: scientific design — complete

- Two scenarios represent two different antibiotic-resistant organisms.
- Prioritize approximate \(R_0\) targets of 1.5 and 2.0 and create the difference through the scenario-specific probability-of-transmission function.
- Use individual realized secondary-case count \(R_i\), as returned in the `rt` column of `get_reproductive_number()`.
- Retain all infected agents, including zero-secondary-case agents.
- Use the provisional incubation, recovery, seeding, and horizon values for the toy prototype.
- Use one compact missingness-augmented R `torch` model that supports both features, \(X\) only, or scenario only.

### Gate 2: calibrated pilot

- Achieved targets within tolerance.
- No invalid transmission probabilities.
- Extraction identities and reproducibility tests pass.
- Runtime and output size projections are acceptable.

### Gate 3: production run

- Freeze configuration and package versions.
- Execute 10,000 simulations per scenario.
- Audit completed replicate IDs and rerun only missing/failed batches.

### Gate 4: ML analysis

- Freeze train/validation/test run IDs before fitting.
- Apply masking only after splitting.
- Report complete-data and each incomplete-data pattern separately.

## 11. Remaining operational question

What compute environment will run the 20,000 production simulations (local workstation, cluster, or cloud), and is there a preferred parallelization system? This does not block construction of the smoke test or calibration pilot.

## 12. Approved working design

- Two separate scenarios representing lower- and higher-transmission antibiotic-resistant organisms.
- Native logit transmission, \(p_{is}=\operatorname{logit}^{-1}(\alpha+\beta X_i+\delta S_s)\), calibrated to produce clearly separated empirical early-epidemic reproduction numbers near 1.5 and 2.0.
- \(R_i\), obtained from `get_reproductive_number()`, as the individual target.
- All infected agents retained, including zeros.
- Redraw \(X\) in each replicate unless later changed.
- Use the provisional epidemiological constants in Section 2 for the smoke test.
- Train with equal probabilities for the three supported observation patterns unless expected real-world feature availability becomes known.
- Compare negative-binomial regression with one compact missingness-augmented R `torch` count model.

## 13. Technical references for the transmission callback

- [`epiworldR` package reference](https://uofuepibio.github.io/epiworldR/reference/index.html)
- [`virus_fun_logit()` and `set_prob_infecting_fun()` documentation](https://uofuepibio.github.io/epiworldR/reference/virus.html)
- [`ModelSEIRCONN` documentation](https://uofuepibio.github.io/epiworldR/reference/ModelSEIRCONN.html)
- [`ModelSEIRCONN` C++ transmission path](https://github.com/UofUEpiBio/epiworld/blob/master/include/epiworld/models/seirconnected.hpp)
- [`virus_fun_logit` C++ implementation](https://github.com/UofUEpiBio/epiworld/blob/master/include/epiworld/virus-meat.hpp)
- [`epiworldR` C++ bindings for virus callbacks](https://github.com/UofUEpiBio/epiworldR/blob/main/src/virus.cpp)
- [R `torch` project and installation guide](https://torch.mlverse.org/)
- [`nn_module()` reference](https://torch.mlverse.org/docs/reference/nn_module.html)
- [`torch_load()` reference](https://torch.mlverse.org/docs/reference/torch_load.html)
