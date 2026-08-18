#' Simulate one heterogeneous SEIR epidemic
#'
#' @param config Named simulation-configuration list.
#' @param scenario Character scenario label.
#' @param replicate_id Positive integer replicate identifier.
#' @param seed Optional integer random seed.
#' @param keep_transmissions Whether to retain the transmission-edge table.
#'
#' @return A list containing agent outcomes, a run summary, and optionally
#'   transmission edges.
#' @export
simulate_one <- function(
    config,
    scenario,
    replicate_id,
    seed               = NULL,
    keep_transmissions = FALSE
) {
  if (is.null(seed)) {
    scenario_offset <- if (scenario == "higher") 100000000L else 0L
    seed <- as.integer(config$base_seed + scenario_offset + replicate_id)
  }
  set.seed(seed)
  x <- stats::runif(config$n_agents)
  bundle <- build_seirconn_model(config, scenario, x)
  epiworldR::run(bundle$model, ndays = config$max_days, seed = seed + 1L)

  result <- extract_simulation_results(
    bundle,
    config = config,
    run_id = sprintf("%s_%06d", scenario, replicate_id),
    seed = seed
  )
  if (!keep_transmissions) result$transmissions <- NULL
  result
}

#' Run an in-memory simulation study
#'
#' @param config Named simulation-configuration list.
#' @param n_reps Number of replicates per scenario.
#' @param scenarios Character vector of scenario labels.
#' @param workers Number of parallel workers.
#' @param keep_transmissions Whether to retain transmission-edge tables.
#'
#' @return A list containing combined agent and run tables and, when requested,
#'   transmission edges.
#' @export
run_simulation_study <- function(
    config,
    n_reps,
    scenarios          = c("lower", "higher"),
    workers            = 1L,
    keep_transmissions = FALSE
) {
  tasks <- expand.grid(
    replicate_id = seq_len(n_reps),
    scenario = scenarios,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  run_task <- function(i) {
    simulate_one(
      config = config,
      scenario = tasks$scenario[i],
      replicate_id = tasks$replicate_id[i],
      keep_transmissions = keep_transmissions
    )
  }

  if (workers > 1L && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(
      seq_len(nrow(tasks)), run_task, mc.cores = workers, mc.preschedule = TRUE
    )
  } else {
    results <- lapply(seq_len(nrow(tasks)), run_task)
  }

  out <- list(
    agents = do.call(rbind, lapply(results, `[[`, "agents")),
    runs = do.call(rbind, lapply(results, `[[`, "runs"))
  )
  if (keep_transmissions) {
    out$transmissions <- do.call(rbind, Map(
      function(result, task_id) {
        data.frame(run_id = task_id, result$transmissions, check.names = FALSE)
      },
      results,
      sprintf("%s_%06d", tasks$scenario, tasks$replicate_id)
    ))
  }
  rownames(out$agents) <- NULL
  rownames(out$runs) <- NULL
  out
}

#' Run and save a batched simulation study
#'
#' @param config Named simulation-configuration list.
#' @param n_reps Number of replicates per scenario.
#' @param batch_size Number of replicates per saved batch.
#' @param workers Number of parallel workers.
#' @param output_dir Directory for ignored generated RDS files.
#'
#' @return A data frame manifest of generated batch files.
#' @export
run_simulation_batches <- function(
    config,
    n_reps     = 10000L,
    batch_size = 100L,
    workers    = 1L,
    output_dir = "data/derived"
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- list()
  for (scenario in c("lower", "higher")) {
    starts <- seq.int(1L, n_reps, by = batch_size)
    for (start in starts) {
      ids <- start:min(start + batch_size - 1L, n_reps)
      if (workers > 1L && .Platform$OS.type != "windows") {
        tasks <- parallel::mclapply(ids, function(id) {
          simulate_one(config, scenario, id, keep_transmissions = FALSE)
        }, mc.cores = workers, mc.preschedule = TRUE)
      } else {
        tasks <- lapply(ids, function(id) {
          simulate_one(config, scenario, id, keep_transmissions = FALSE)
        })
      }
      batch <- list(
        agents = do.call(rbind, lapply(tasks, `[[`, "agents")),
        runs = do.call(rbind, lapply(tasks, `[[`, "runs")),
        config = config
      )
      path <- file.path(
        output_dir,
        sprintf("simulation_%s_%06d_%06d.rds", scenario, min(ids), max(ids))
      )
      saveRDS(batch, path, compress = "xz")
      manifest[[length(manifest) + 1L]] <- data.frame(
        scenario = scenario,
        first_replicate = min(ids),
        last_replicate = max(ids),
        path = path,
        stringsAsFactors = FALSE
      )
    }
  }
  manifest <- do.call(rbind, manifest)
  saveRDS(manifest, file.path(output_dir, "manifest.rds"))
  manifest
}
