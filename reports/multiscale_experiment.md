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
design_summary <- aggregate(
  profile_id ~ site_id,
  study$design,
  length
)
names(design_summary)[2] <- "profiles"
knitr::kable(design_summary, row.names = FALSE)
```

| site_id | profiles |
|:--------|---------:|
| site_a  |       25 |
| site_b  |       25 |
| site_c  |       25 |
| site_d  |       25 |

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
| site_a_001 | site_a | 14 | 0 | 7.604 | FALSE | -0.037 | -0.495 | TRUE | 1 | 1 | 1 | 1 |
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
parameter_summary <- do.call(rbind, lapply(
  split(parameter_recovery, parameter_recovery$parameter),
  function(rows) data.frame(
    parameter = rows$parameter[1],
    profiles = nrow(rows),
    mean_truth = mean(rows$truth),
    mean_fitted = mean(rows$fitted),
    rmse = sqrt(mean((rows$fitted - rows$truth)^2)),
    boundary_fraction = mean(rows$boundary)
  )
))
knitr::kable(parameter_summary, digits = 3, row.names = FALSE)
```

| parameter | profiles | mean_truth | mean_fitted |  rmse | boundary_fraction |
|:----------|---------:|-----------:|------------:|------:|------------------:|
| delta     |      100 |      0.134 |       0.125 | 0.066 |              0.92 |
| gamma     |      100 |      0.086 |       0.095 | 0.057 |              0.92 |
| h         |      100 |      0.149 |       0.161 | 0.088 |              0.92 |

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
| bg1 | 0.000 | 0.570 | 0.206 | 0.217 |
| bg2 | 0.464 | 0.000 | 0.175 | 0.185 |
| bg3 | 0.123 | 0.128 | 0.000 | 0.398 |
| bg4 | 0.098 | 0.102 | 0.301 | 0.000 |

Normalized transfer exposure

``` r
knitr::kable(example$re$matrix, digits = 3, caption = "Effective within-host next-generation matrix")
```

|     |   bg1 |   bg2 |   bg3 |   bg4 |
|:----|------:|------:|------:|------:|
| bg1 | 0.000 | 0.329 | 0.134 | 0.141 |
| bg2 | 0.260 | 0.000 | 0.108 | 0.113 |
| bg3 | 0.071 | 0.074 | 0.000 | 0.212 |
| bg4 | 0.052 | 0.055 | 0.152 | 0.000 |

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
comparison_summary <- aggregate(
  relative_overstatement ~ site_id,
  comparison,
  function(value) c(mean = mean(value), maximum = max(value))
)
comparison_summary <- data.frame(
  site_id = comparison_summary$site_id,
  mean_overstatement = comparison_summary$relative_overstatement[, "mean"],
  maximum_overstatement = comparison_summary$relative_overstatement[, "maximum"]
)
knitr::kable(comparison_summary, digits = 3, row.names = FALSE)
```

| site_id | mean_overstatement | maximum_overstatement |
|:--------|-------------------:|----------------------:|
| site_a  |              0.016 |                 0.074 |
| site_b  |              0.044 |                 0.116 |
| site_c  |              0.051 |                 0.112 |
| site_d  |              0.035 |                 0.085 |

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
|               281 |                   276 |                      0.982 |

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
| baseline               |     0.704 |     1.059 |      5.222 |      5.413 |
| conjugation_inhibition |     0.410 |     0.688 |      2.135 |      2.491 |
| shorter_antibiotic     |     0.704 |     0.859 |      5.222 |      4.929 |

## Missing-modality emulator

The emulator is trained on fitted mechanistic targets from 100 synthetic
host profiles and three intervention settings per profile. During model
development, one model is evaluated on an entirely unseen synthetic
site. The browser model is then retrained with every site represented
and evaluated on profiles held out within each site, because the
deployed interface supports all four site profiles. Errors against
fitted targets measure emulation; differences between fitted and truth
targets measure mechanistic fitting error. Neither is a Bayesian
uncertainty interval or a confidence interval.

``` r
evaluation <- readRDS("../artifacts/multiscale_emulator_evaluation.rds")
complete_input <- evaluation[evaluation$pattern == "all", ]
knitr::kable(complete_input, digits = 3, row.names = FALSE)
```

| evaluation_scope | pattern | target     |   n |  rmse | balanced_accuracy |
|:-----------------|:--------|:-----------|----:|------:|------------------:|
| leave_site_out   | all     | r0_within  |  75 | 0.827 |             0.500 |
| leave_site_out   | all     | re_within  |  75 | 1.124 |             0.500 |
| leave_site_out   | all     | r0_between |  75 | 7.034 |             0.500 |
| leave_site_out   | all     | re_between |  75 | 6.467 |             0.500 |
| profile_holdout  | all     | r0_within  |  36 | 0.272 |             0.667 |
| profile_holdout  | all     | re_within  |  36 | 0.353 |             0.803 |
| profile_holdout  | all     | r0_between |  36 | 2.233 |             0.976 |
| profile_holdout  | all     | re_between |  36 | 1.856 |             0.978 |

``` r
profile_missingness <- evaluation[
  evaluation$evaluation_scope == "profile_holdout" &
    evaluation$target %in% c("re_within", "re_between"),
]
knitr::kable(profile_missingness, digits = 3, row.names = FALSE)
```

| evaluation_scope | pattern           | target     |   n |  rmse | balanced_accuracy |
|:-----------------|:------------------|:-----------|----:|------:|------------------:|
| profile_holdout  | all               | re_within  |  36 | 0.353 |             0.803 |
| profile_holdout  | all               | re_between |  36 | 1.856 |             0.978 |
| profile_holdout  | no_quantitative   | re_within  |  36 | 0.510 |             0.709 |
| profile_holdout  | no_quantitative   | re_between |  36 | 3.322 |             0.875 |
| profile_holdout  | no_genomic        | re_within  |  36 | 0.595 |             0.666 |
| profile_holdout  | no_genomic        | re_between |  36 | 3.547 |             0.747 |
| profile_holdout  | no_clinical       | re_within  |  36 | 0.542 |             0.615 |
| profile_holdout  | no_clinical       | re_between |  36 | 3.717 |             0.846 |
| profile_holdout  | quantitative_only | re_within  |  36 | 0.626 |             0.550 |
| profile_holdout  | quantitative_only | re_between |  36 | 4.079 |             0.615 |
| profile_holdout  | genomic_only      | re_within  |  36 | 0.673 |             0.435 |
| profile_holdout  | genomic_only      | re_between |  36 | 4.774 |             0.846 |
| profile_holdout  | clinical_only     | re_within  |  36 | 0.709 |             0.550 |
| profile_holdout  | clinical_only     | re_between |  36 | 4.528 |             0.600 |

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
| r0_within  | 0.368 |
| re_within  | 0.506 |
| r0_between | 3.994 |
| re_between | 3.409 |

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
- Leave-site-out performance remains weak because the synthetic sites
  differ in both biological parameter distributions and structural
  modality availability; the browser model is a functional
  demonstration, not a validated predictor for an unseen setting.
- The four-versus-five Peru sample schedule, Hi-C detection limit, and
  transfer subscript inconsistency require confirmation before empirical
  use.
