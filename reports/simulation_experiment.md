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
| Infected while at least 95% of the population was susceptible | `early_susceptible_fraction` | Rules out counts capped by susceptible depletion |
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
| higher | 9613099 | 3654073 | 9144999 | 35531 | 467186 | 0.049 |
| lower | 852123 | 740604 | 485557 | 240207 | 111471 | 0.131 |

``` r
analysis_agents <- filter_analysis_agents(study$agents)
```

The two scenarios lose agents to different criteria. The `higher`
epidemic leaves its early phase after about five weeks and then infects
essentially the whole population, so susceptible depletion is what
excludes most of its agents. The `lower` epidemic stays inside the early
phase far longer, and the day-60 infection cutoff is what excludes most
of its agents. Both scenarios are supercritical, so both contribute a
well-populated analysis window: eligible `higher` agents outnumber
eligible `lower` agents by roughly four to one rather than by the two
orders of magnitude a subcritical `lower` scenario would produce.

Most runs still have exposed or infectious agents at the horizon. That
is a statement about the epidemic, not about the analyzed agents: an
agent enters the window only when its own infectious period has ended,
and an agent infected on or before day 60 has at least 60 days to
recover before day 120.

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
| higher   |                    34 |              9625 |                    999 |
| lower    |                   108 |               720 |                    859 |

## Why the unfiltered average is bounded by one

Every non-seed infection is counted exactly once as somebody’s secondary
case. Averaging `R_i` over all infected agents therefore returns
`(infections - seeds) / infections`, which sits just below one
regardless of how transmissible the organism is. The analysis window
removes that artifact.

``` r
comparison <- do.call(rbind, lapply(
  split(study$agents, study$agents$scenario),
  function(data) {
    within_cutoff <- data[
      data$infection_day <= config$infection_cutoff_day,
    ]
    data.frame(
      scenario = data$scenario[[1]],
      all_infected = mean(data$secondary_cases),
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
| higher   |        0.999 |              1.362 |             2.987 |       3.0 |
| lower    |        0.994 |              1.477 |             1.477 |       1.5 |

The first column is the identity. The second still includes agents
infected after the susceptible pool began to deplete, which is why it
stays far below the target in the `higher` scenario; in the `lower`
scenario the early phase outlasts the day-60 cutoff, so the second and
third columns coincide. Only the third column estimates the
transmissibility the calibration targeted.

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
| higher   |          467186 |   2.987 | 3.453 |   0 |   1 |      2 |   4 |  10 |         0.218 |
| lower    |          111471 |   1.477 | 2.015 |   0 |   0 |      1 |   2 |   5 |         0.399 |

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
| higher   |   2.987 | 0.005 |          467186 |
| lower    |   1.477 | 0.006 |          111471 |

``` r
knitr::kable(run_early, digits = 3)
```

| scenario | mean_run_ri | sd_run_ri | runs |
|:---------|------------:|----------:|-----:|
| higher   |       2.986 |     0.166 | 1000 |
| lower    |       1.301 |     0.383 | 1000 |

The pooled means sit close to the 1.5 and 3 targets. The run-level
`lower` mean is smaller than the pooled mean because every run counts
once no matter how many eligible agents it contributed, and the number
of eligible agents a run produces is itself correlated with that run’s
realized transmission. These values vary around the targets because this
is a finite stochastic model. The calibration script refines the
scenario intercepts before the production run.
