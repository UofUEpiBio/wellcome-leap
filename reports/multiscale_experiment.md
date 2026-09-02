# AMPLIFY multiscale pARG experiment


This synthetic experiment connects longitudinal omic observations to a
four-background within-host ODE, explicit reference and effective
reproduction numbers, a summary-coupled carriage ABM, and a
missing-modality emulator. It models pARG carriage and acquisition, not
clinical infection.

## Generated study

``` r
source("../config/multiscale.R")
source("../R/within_host.R")
source("../R/omics.R")
source("../R/reproduction.R")
source("../R/multiscale.R")

study <- readRDS("../data/derived/multiscale/fitted_study.rds")
truth_targets <- combine_multiscale_targets(study$truth)
fitted_targets <- combine_multiscale_targets(study$fitted)
knitr::kable(study$design, row.names = FALSE)
```

| profile_number | site_id | profile_id |
|---------------:|:--------|:-----------|
|              1 | site_a  | site_a_001 |
|              2 | site_a  | site_a_002 |
|              1 | site_b  | site_b_001 |
|              2 | site_b  | site_b_002 |
|              1 | site_c  | site_c_001 |
|              2 | site_c  | site_c_002 |
|              1 | site_d  | site_d_001 |
|              2 | site_d  | site_d_002 |

Each of four sites changes the distribution of the three fitted
parameters, contact rate, susceptible fraction, and structural assay
availability. Every host is observed on days 0, 3, 7, 14, and 30.

``` r
knitr::kable(study$simulated[[1]]$observations, digits = 3, row.names = FALSE)
```

| profile_id | site_id | day | antibiotic_active | qpcr_log10 | qpcr_censored | ecoli_log10 | klebsiella_log10 | culture_positive | link_bg1 | link_bg2 | link_bg3 | link_bg4 |
|:---|:---|---:|---:|---:|:---|---:|---:|:---|---:|---:|---:|---:|
| site_a_001 | site_a | 0 | 1 | 7.682 | FALSE | -0.301 | -0.670 | TRUE | 1 | 0 | 0 | 0 |
| site_a_001 | site_a | 3 | 1 | 7.248 | FALSE | -0.271 | -0.287 | TRUE | 1 | 0 | 0 | 1 |
| site_a_001 | site_a | 7 | 0 | 7.623 | FALSE | -0.200 | -0.475 | TRUE | 1 | 1 | 1 | 1 |
| site_a_001 | site_a | 14 | 0 | 7.605 | FALSE | -0.037 | -0.495 | TRUE | 1 | 1 | 1 | 1 |
| site_a_001 | site_a | 30 | 0 | 7.409 | FALSE | -0.238 | -0.393 | TRUE | 1 | 1 | 1 | 0 |

## Truth, observations, and parameter fitting

The ODE uses two *E. coli* and two *Klebsiella* genomic backgrounds. All
four backgrounds share conjugation `h`, division-associated segregation
`gamma`, and antibiotic selection `delta`; a fixed matrix supplies
relative within- and cross-species conjugation rates.

``` r
parameter_recovery <- do.call(rbind, Map(
  function(truth, fitted) {
    data.frame(
      profile_id = truth$profile$profile_id,
      site_id = truth$profile$site_id,
      parameter = names(truth$parameters),
      truth = unname(truth$parameters),
      fitted = unname(fitted$parameters),
      converged = fitted$fit$convergence == 0,
      boundary = fitted$fit$at_boundary
    )
  },
  study$truth,
  study$fitted
))
parameter_recovery$absolute_error <- abs(
  parameter_recovery$fitted - parameter_recovery$truth
)
knitr::kable(parameter_recovery, digits = 3, row.names = FALSE)
```

| profile_id | site_id | parameter | truth | fitted | converged | boundary | absolute_error |
|:-----------|:--------|:----------|------:|-------:|:----------|:---------|---------------:|
| site_a_001 | site_a  | h         | 0.104 |  0.046 | TRUE      | FALSE    |          0.058 |
| site_a_001 | site_a  | gamma     | 0.124 |  0.060 | TRUE      | FALSE    |          0.064 |
| site_a_001 | site_a  | delta     | 0.096 |  0.070 | TRUE      | FALSE    |          0.026 |
| site_a_002 | site_a  | h         | 0.162 |  0.252 | TRUE      | TRUE     |          0.091 |
| site_a_002 | site_a  | gamma     | 0.082 |  0.160 | TRUE      | TRUE     |          0.078 |
| site_a_002 | site_a  | delta     | 0.121 |  0.178 | TRUE      | TRUE     |          0.058 |
| site_b_001 | site_b  | h         | 0.091 |  0.041 | TRUE      | FALSE    |          0.051 |
| site_b_001 | site_b  | gamma     | 0.109 |  0.061 | TRUE      | FALSE    |          0.048 |
| site_b_001 | site_b  | delta     | 0.154 |  0.161 | TRUE      | FALSE    |          0.007 |
| site_b_002 | site_b  | h         | 0.159 |  0.090 | TRUE      | FALSE    |          0.068 |
| site_b_002 | site_b  | gamma     | 0.089 |  0.045 | TRUE      | FALSE    |          0.045 |
| site_b_002 | site_b  | delta     | 0.105 |  0.163 | TRUE      | FALSE    |          0.059 |
| site_c_001 | site_c  | h         | 0.140 |  0.171 | TRUE      | TRUE     |          0.031 |
| site_c_001 | site_c  | gamma     | 0.086 |  0.160 | TRUE      | TRUE     |          0.074 |
| site_c_001 | site_c  | delta     | 0.138 |  0.156 | TRUE      | TRUE     |          0.017 |
| site_c_002 | site_c  | h         | 0.119 |  0.117 | TRUE      | TRUE     |          0.002 |
| site_c_002 | site_c  | gamma     | 0.056 |  0.040 | TRUE      | TRUE     |          0.016 |
| site_c_002 | site_c  | delta     | 0.092 |  0.057 | TRUE      | TRUE     |          0.035 |
| site_d_001 | site_d  | h         | 0.137 |  0.350 | TRUE      | TRUE     |          0.213 |
| site_d_001 | site_d  | gamma     | 0.074 |  0.160 | TRUE      | TRUE     |          0.086 |
| site_d_001 | site_d  | delta     | 0.203 |  0.040 | TRUE      | TRUE     |          0.163 |
| site_d_002 | site_d  | h         | 0.155 |  0.033 | TRUE      | TRUE     |          0.121 |
| site_d_002 | site_d  | gamma     | 0.075 |  0.040 | TRUE      | TRUE     |          0.035 |
| site_d_002 | site_d  | delta     | 0.163 |  0.250 | TRUE      | TRUE     |          0.087 |

``` r
groups <- split(parameter_recovery, parameter_recovery$parameter)
matplot(
  do.call(cbind, lapply(groups, `[[`, "truth")),
  do.call(cbind, lapply(groups, `[[`, "fitted")),
  pch = 19,
  col = c("#1b9e77", "#d95f02", "#7570b3"),
  xlab = "Truth",
  ylab = "Fitted value"
)
abline(0, 1, lty = 2)
legend(
  "topleft",
  legend = names(groups),
  col = c("#1b9e77", "#d95f02", "#7570b3"),
  pch = 19,
  bty = "n"
)
```

![Fitted versus true shared within-host
parameters.](multiscale_experiment_files/figure-commonmark/unnamed-chunk-4-1.png)

## Within-host reproduction

For donor background `j`, the implementation evolves only the introduced
donor lineage around the pARG-free susceptible trajectory. Recipient
acquisitions do not become donors inside that calculation. The
normalized exposure and first-generation matrix are

$$
\Lambda_{ij}=q^{-1}\int h\omega_{ij}S_i^0(t)R_j^{(1)}(t)/K\,dt,
\qquad
K_{ij}=1-\exp[-p_{est,i}\Lambda_{ij}].
$$

The seed normalization makes the result invariant to the arbitrary
numerical inoculum. Diagonal entries are zero because an already
positive background is not a new genomic background.

``` r
example <- study$fitted[[1]]$baseline$within
knitr::kable(example$re$lambda, digits = 3, caption = "Normalized transfer exposure")
```

|     |   bg1 |   bg2 |   bg3 |   bg4 |
|:----|------:|------:|------:|------:|
| bg1 | 0.000 | 0.571 | 0.206 | 0.217 |
| bg2 | 0.465 | 0.000 | 0.175 | 0.185 |
| bg3 | 0.123 | 0.128 | 0.000 | 0.398 |
| bg4 | 0.098 | 0.102 | 0.301 | 0.000 |

Normalized transfer exposure

``` r
knitr::kable(example$re$matrix, digits = 3, caption = "Effective within-host next-generation matrix")
```

|     |   bg1 |   bg2 |   bg3 |   bg4 |
|:----|------:|------:|------:|------:|
| bg1 | 0.000 | 0.329 | 0.134 | 0.141 |
| bg2 | 0.261 | 0.000 | 0.108 | 0.113 |
| bg3 | 0.071 | 0.074 | 0.000 | 0.212 |
| bg4 | 0.052 | 0.055 | 0.153 | 0.000 |

Effective within-host next-generation matrix

``` r
data.frame(
  metric = c("R0 within host", "Re within host"),
  value = c(example$r0$reproduction, example$re$reproduction)
) |>
  knitr::kable(digits = 3, row.names = FALSE)
```

| metric         | value |
|:---------------|------:|
| R0 within host | 0.303 |
| Re within host | 0.421 |

## Between-host reproduction

The between-host expectation integrates the complete resistant-fraction
trajectory rather than applying the inoculum function to its mean:

$$
R_{e,i}^{BH}=\int_0^{D_i}c_i s_i p_{est,i}
\{1-\exp[-\kappa p_i(t)]\}\,dt.
$$

The product-form comparator is slightly larger, as expected from the
concavity of the inoculum function.

``` r
comparison <- fitted_targets[
  fitted_targets$intervention == "baseline",
  c("profile_id", "site_id", "re_between", "product_between")
]
comparison$relative_overstatement <-
  comparison$product_between / comparison$re_between - 1
knitr::kable(comparison, digits = 3, row.names = FALSE)
```

| profile_id | site_id | re_between | product_between | relative_overstatement |
|:-----------|:--------|-----------:|----------------:|-----------------------:|
| site_a_001 | site_a  |      0.567 |           0.569 |                  0.004 |
| site_a_002 | site_a  |      4.989 |           5.131 |                  0.029 |
| site_b_001 | site_b  |      0.949 |           0.959 |                  0.010 |
| site_b_002 | site_b  |      5.587 |           5.834 |                  0.044 |
| site_c_001 | site_c  |      1.713 |           1.715 |                  0.002 |
| site_c_002 | site_c  |      8.110 |           8.964 |                  0.105 |
| site_d_001 | site_d  |     10.613 |          11.199 |                  0.055 |
| site_d_002 | site_d  |      3.153 |           3.173 |                  0.006 |

The stochastic ABM retains realized acquisitions as a diagnostic rather
than calling one offspring count a reproduction number.

``` r
abm_summary <- data.frame(
  infected_carriers = nrow(study$abm$agents),
  realized_acquisitions = sum(study$abm$agents$realized_secondary_acquisitions),
  mean_realized_acquisitions = mean(
    study$abm$agents$realized_secondary_acquisitions
  )
)
knitr::kable(abm_summary, digits = 3, row.names = FALSE)
```

| infected_carriers | realized_acquisitions | mean_realized_acquisitions |
|------------------:|----------------------:|---------------------------:|
|               104 |                    99 |                      0.952 |

## Paired intervention scenarios

Shorter antibiotic exposure changes the current-condition values but not
the antibiotic-free reference. Conjugation inhibition multiplies `h` by
0.5 and therefore changes both reference and current-condition values.

``` r
intervention_summary <- aggregate(
  cbind(r0_within, re_within, r0_between, re_between) ~ intervention,
  fitted_targets,
  mean
)
knitr::kable(intervention_summary, digits = 3, row.names = FALSE)
```

| intervention           | r0_within | re_within | r0_between | re_between |
|:-----------------------|----------:|----------:|-----------:|-----------:|
| baseline               |     0.540 |     0.896 |      3.911 |      4.460 |
| conjugation_inhibition |     0.294 |     0.519 |      0.676 |      1.157 |
| shorter_antibiotic     |     0.540 |     0.678 |      3.911 |      3.828 |

## Missing-modality emulator

The emulator is trained on fitted mechanistic targets and tested on an
unseen synthetic site. Errors against fitted targets measure emulation;
differences between fitted and truth targets measure mechanistic fitting
error. Neither is a Bayesian uncertainty interval.

``` r
evaluation <- readRDS("../artifacts/multiscale_emulator_evaluation.rds")
knitr::kable(evaluation, digits = 3, row.names = FALSE)
```

| pattern           | target     |   n |  rmse | balanced_accuracy |
|:------------------|:-----------|----:|------:|------------------:|
| all               | r0_within  |   6 | 0.366 |             1.000 |
| all               | re_within  |   6 | 0.278 |             0.500 |
| all               | r0_between |   6 | 5.605 |             0.833 |
| all               | re_between |   6 | 2.832 |             1.000 |
| no_quantitative   | r0_within  |   6 | 0.366 |             1.000 |
| no_quantitative   | re_within  |   6 | 0.278 |             0.500 |
| no_quantitative   | r0_between |   6 | 5.605 |             0.833 |
| no_quantitative   | re_between |   6 | 2.832 |             1.000 |
| no_genomic        | r0_within  |   6 | 0.363 |             1.000 |
| no_genomic        | re_within  |   6 | 0.248 |             0.500 |
| no_genomic        | r0_between |   6 | 6.207 |             0.500 |
| no_genomic        | re_between |   6 | 4.228 |             1.000 |
| no_clinical       | r0_within  |   6 | 0.176 |             1.000 |
| no_clinical       | re_within  |   6 | 0.227 |             0.500 |
| no_clinical       | r0_between |   6 | 6.186 |             1.000 |
| no_clinical       | re_between |   6 | 5.087 |             0.500 |
| quantitative_only | r0_within  |   6 | 0.265 |             1.000 |
| quantitative_only | re_within  |   6 | 0.207 |             0.500 |
| quantitative_only | r0_between |   6 | 6.539 |             0.500 |
| quantitative_only | re_between |   6 | 5.492 |             1.000 |
| genomic_only      | r0_within  |   6 | 0.176 |             1.000 |
| genomic_only      | re_within  |   6 | 0.227 |             0.500 |
| genomic_only      | r0_between |   6 | 6.186 |             1.000 |
| genomic_only      | re_between |   6 | 5.087 |             0.500 |
| clinical_only     | r0_within  |   6 | 0.363 |             1.000 |
| clinical_only     | re_within  |   6 | 0.248 |             0.500 |
| clinical_only     | r0_between |   6 | 6.207 |             0.500 |
| clinical_only     | re_between |   6 | 4.228 |             1.000 |

``` r
matched <- merge(
  truth_targets,
  fitted_targets,
  by = c("profile_id", "site_id", "intervention"),
  suffixes = c("_truth", "_fitted")
)
targets <- c("r0_within", "re_within", "r0_between", "re_between")
fitting_error <- data.frame(
  target = targets,
  rmse = vapply(targets, function(target) {
    sqrt(mean((
      matched[[paste0(target, "_fitted")]] -
        matched[[paste0(target, "_truth")]]
    )^2))
  }, numeric(1))
)
knitr::kable(fitting_error, digits = 3, row.names = FALSE)
```

| target     |  rmse |
|:-----------|------:|
| r0_within  | 0.163 |
| re_within  | 0.240 |
| r0_between | 3.049 |
| re_between | 2.395 |

## Scientific limitations exposed by the prototype

- The proposal must define and normalize the within-host infectious
  unit.
- Genomic host switches are interval-censored evidence, not direct
  observations of conjugation.
- Five observations do not identify background-specific biological
  parameters; this prototype therefore fits only three shared values.
- Carriage and clinical infection are different outcomes.
- Deterministic ODEs do not represent rare-transfer extinction.
- Point fitting does not propagate mechanistic uncertainty into the
  emulator.
- The four-versus-five Peru sample schedule, Hi-C detection limit, and
  transfer subscript inconsistency require confirmation before empirical
  use.
