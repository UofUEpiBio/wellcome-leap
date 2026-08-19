# AMPLIFY simulation experiment


This report summarizes the full production study for the two
antibiotic-resistant organism scenarios. It uses 1,000 independent
simulations per scenario with 10,000 agents per simulation and a 120-day
horizon. Each simulation draws the agent feature from a standard normal
distribution.

## Setup and execution

``` r
source("../config/simulation.R")
source("../R/build_model.R")
source("../R/extract_outcomes.R")
source("../R/simulate.R")

config <- default_simulation_config()
study <- load_simulation_batches(
  manifest_path = "../data/derived/production/manifest.rds",
  root_dir       = ".."
)
run_counts <- table(study$runs$scenario)
if (!identical(as.integer(run_counts[c("lower", "higher")]), c(1000L, 1000L))) {
  stop("The production manifest must contain 1,000 runs per scenario.")
}
knitr::kable(data.frame(scenario = names(run_counts), runs = as.integer(run_counts)))
```

| scenario | runs |
|:---------|-----:|
| higher   | 1000 |
| lower    | 1000 |

## Analysis window

A realized individual reproduction number is only informative about the
agent who produced it when that agent had a full opportunity to transmit
and was infected while the susceptible pool was still essentially
intact. Four criteria implement this, and an agent must satisfy all of
them to be `analysis_eligible`:

| Criterion | Configuration | Reason |
|----|----|----|
| Infected on or before day 60 | `infection_cutoff_day` | Leaves at least 60 days to transmit inside the 120-day horizon |
| Transmission window of at least 30 days | `min_transmission_days` | Rules out right-censored offspring counts |
| Infected while at least 90% of the population was susceptible | `early_susceptible_fraction` | Rules out counts capped by susceptible depletion |
| Not exposed or infectious at the horizon | `max_days` | Rules out agents whose infectious period never ended |

The pseudo-source row identified by `source = -1` is excluded before
individual outcomes are constructed. Eligible agents include those with
zero secondary infections.

``` r
eligibility <- do.call(rbind, lapply(
  split(study$agents, study$agents$scenario),
  function(data) {
    data.frame(
      scenario = data$scenario[[1]],
      infected_agents = nrow(data),
      infected_after_cutoff = sum(
        data$infection_day > config$infection_cutoff_day
      ),
      infected_after_early_phase = sum(!data$early_phase),
      censored_at_horizon = sum(!data$outcome_complete),
      eligible_agents = sum(data$analysis_eligible),
      eligible_fraction = mean(data$analysis_eligible)
    )
  }
))
rownames(eligibility) <- NULL
knitr::kable(eligibility, digits = 3)
```

| scenario | infected_agents | infected_after_cutoff | infected_after_early_phase | censored_at_horizon | eligible_agents | eligible_fraction |
|:---|---:|---:|---:|---:|---:|---:|
| higher | 9992710 | 10734 | 9100964 | 89 | 891745 | 0.089 |
| lower | 9128 | 153 | 0 | 7 | 8975 | 0.983 |

``` r
analysis_agents <- filter_analysis_agents(study$agents)
```

Susceptible depletion, not the infection-day cutoff, is what excludes
most `higher`-scenario agents. That scenario leaves its early phase
after roughly three weeks and then goes on to infect essentially the
whole population well before day 60, so nearly all of its infections
happen after susceptible depletion has begun.

``` r
run_shape <- do.call(rbind, lapply(
  split(study$runs, study$runs$scenario),
  function(data) {
    data.frame(
      scenario = data$scenario[[1]],
      median_early_last_day = stats::median(data$early_last_day),
      median_final_size = stats::median(data$final_epidemic_size),
      runs_active_at_horizon = sum(!data$outcome_complete)
    )
  }
))
rownames(run_shape) <- NULL
knitr::kable(run_shape, digits = 1)
```

| scenario | median_early_last_day | median_final_size | runs_active_at_horizon |
|:---------|----------------------:|------------------:|-----------------------:|
| higher   |                    22 |              9993 |                     76 |
| lower    |                   120 |                 7 |                      4 |

## Why the unfiltered average is bounded by one

Every non-seed infection is counted exactly once as somebody’s secondary
case. Averaging `R_i` over all infected agents therefore returns
`(final size - seeds) / final size`, which approaches one from below in
any completed epidemic regardless of how transmissible the organism is.
The analysis window removes that artifact.

``` r
comparison <- do.call(rbind, lapply(
  split(study$agents, study$agents$scenario),
  function(data) {
    complete <- data[data$outcome_complete, ]
    within_cutoff <- complete[
      complete$infection_day <= config$infection_cutoff_day,
    ]
    data.frame(
      scenario = data$scenario[[1]],
      all_infected = mean(complete$secondary_cases),
      infected_by_cutoff = mean(within_cutoff$secondary_cases),
      analysis_eligible = mean(
        data$secondary_cases[data$analysis_eligible]
      ),
      target_r0 = config$target_r0[[data$scenario[[1]]]]
    )
  }
))
rownames(comparison) <- NULL
knitr::kable(comparison, digits = 3)
```

| scenario | all_infected | infected_by_cutoff | analysis_eligible | target_r0 |
|:---------|-------------:|-------------------:|------------------:|----------:|
| higher   |        1.000 |              1.000 |             4.012 |       4.0 |
| lower    |        0.452 |              0.453 |             0.453 |       0.5 |

## Individual reproduction-number distribution

The target is each eligible agent’s realized number of secondary
infections, `R_i`, including zero.

``` r
#' Summarize an individual reproduction-number sample
#'
#' @param data Agent-level simulation data for one scenario.
#'
#' @return A one-row data frame of distribution summaries.
summarize_ri <- function(data) {
  quantiles <- stats::quantile(
    data$secondary_cases,
    probs = c(0.05, 0.25, 0.50, 0.75, 0.95)
  )
  data.frame(
    infected_agents = nrow(data),
    mean_ri = mean(data$secondary_cases),
    sd_ri = stats::sd(data$secondary_cases),
    p05 = quantiles[[1]],
    p25 = quantiles[[2]],
    median = quantiles[[3]],
    p75 = quantiles[[4]],
    p95 = quantiles[[5]],
    zero_fraction = mean(data$secondary_cases == 0)
  )
}

ri_summary <- do.call(
  rbind,
  lapply(split(analysis_agents, analysis_agents$scenario), summarize_ri)
)
ri_summary$scenario <- rownames(ri_summary)
rownames(ri_summary) <- NULL
ri_summary <- ri_summary[c("scenario", setdiff(names(ri_summary), "scenario"))]
knitr::kable(ri_summary, digits = 3)
```

| scenario | infected_agents | mean_ri | sd_ri | p05 | p25 | median | p75 | p95 | zero_fraction |
|:---------|----------------:|--------:|------:|----:|----:|-------:|----:|----:|--------------:|
| higher   |          891745 |   4.012 | 5.066 |   0 |   1 |      2 |   6 |  14 |         0.233 |
| lower    |            8975 |   0.453 | 1.061 |   0 |   0 |      0 |   1 |   2 |         0.739 |

``` r
boxplot(
  secondary_cases ~ scenario,
  data = analysis_agents,
  outline = FALSE,
  xlab = "Scenario",
  ylab = "Individual secondary infections (Ri)",
  col = c("#6baed6", "#fd8d3c")
)
```

![Distribution of individual secondary infections by organism
scenario.](simulation_experiment_files/figure-commonmark/unnamed-chunk-6-1.png)

The next plot retains the large mass at zero and groups the upper tail
for readability.

``` r
ri_group <- pmin(analysis_agents$secondary_cases, 7L)
distribution <- prop.table(
  table(analysis_agents$scenario, factor(ri_group, levels = 0:7)),
  margin = 1
)
barplot(
  distribution,
  beside = TRUE,
  names.arg = c(as.character(0:6), "7+"),
  legend.text = rownames(distribution),
  args.legend = list(x = "topright", bty = "n"),
  xlab = "Individual secondary infections (Ri)",
  ylab = "Fraction of infected agents",
  col = c("#fd8d3c", "#6baed6")
)
```

![Probability mass of individual secondary infections, with values of 7
or more
grouped.](simulation_experiment_files/figure-commonmark/unnamed-chunk-7-1.png)

## Achieved early-epidemic reproduction

Calibration targets the mean `R_i` of analysis-eligible agents. The
table reports both the pooled eligible-agent mean and the distribution
of run-level eligible means.

``` r
pooled_early <- do.call(rbind, lapply(
  split(analysis_agents, analysis_agents$scenario),
  function(x) {
    data.frame(
      scenario = x$scenario[[1]],
      mean_ri = mean(x$secondary_cases),
      se = stats::sd(x$secondary_cases) / sqrt(nrow(x)),
      infected_agents = nrow(x)
    )
  }
))
run_early <- do.call(rbind, lapply(split(study$runs, study$runs$scenario), function(x) {
  data.frame(
    scenario = x$scenario[[1]],
    mean_run_ri = mean(x$eligible_mean_ri, na.rm = TRUE),
    sd_run_ri = stats::sd(x$eligible_mean_ri, na.rm = TRUE),
    runs = nrow(x)
  )
}))
rownames(pooled_early) <- NULL
rownames(run_early) <- NULL
knitr::kable(pooled_early, digits = 3)
```

| scenario | mean_ri |    se | infected_agents |
|:---------|--------:|------:|----------------:|
| higher   |   4.012 | 0.005 |          891745 |
| lower    |   0.453 | 0.011 |            8975 |

``` r
knitr::kable(run_early, digits = 3)
```

| scenario | mean_run_ri | sd_run_ri | runs |
|:---------|------------:|----------:|-----:|
| higher   |       4.018 |     0.136 | 1000 |
| lower    |       0.328 |     0.253 | 1000 |

The pooled means sit close to the 0.5 and 4.0 targets. The run-level
`lower` mean is smaller than the pooled mean because most `lower` runs
go extinct with only their five seeds, so each of those runs contributes
a mean of zero with equal weight. These values vary around the targets
because this is a finite stochastic model. The calibration script
refines the scenario intercepts before the production run.
