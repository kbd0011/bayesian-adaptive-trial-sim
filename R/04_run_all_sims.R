# -----------------------------------------------------------------------------
# Purpose : Orchestrate all trial simulations across 6 scenarios x 2 designs.
#           Dispatches each (scenario, design, sim_id) to the appropriate
#           simulator via furrr; saves raw per-trial results to one RDS per
#           (scenario, design) pair.
# Inputs  : config.yml + outputs/sims/design_params.rds
# Outputs : outputs/sims/raw_<scenario>_<design>.rds (one per pair, 12 total)
# CLI     : --n-sims N       override CONFIG$simulation$n_sims_full
#           --workers N      override CONFIG$simulation$parallel_workers
#           --shard-id K     1-based shard index for matrix sharding
#           --total-shards N total number of shards (default 1 = no sharding)
# Usage   : Rscript R/04_run_all_sims.R --n-sims 100
#           Rscript R/04_run_all_sims.R --n-sims 10000 --shard-id 3 --total-shards 10
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))
source(here::here("R/02_sim_fixed.R"))
source(here::here("R/03_sim_adaptive.R"))

# --- parse CLI args ----------------------------------------------------------
.parse_arg <- function(args, key, default = NULL, as = as.integer) {
  i <- which(args == key)
  if (length(i) == 1L && i < length(args)) as(args[i + 1L]) else default
}
cli_args     <- commandArgs(trailingOnly = TRUE)
n_sims       <- .parse_arg(cli_args, "--n-sims",
                           default = CONFIG$simulation$n_sims_full)
n_workers    <- .parse_arg(cli_args, "--workers",
                           default = CONFIG$simulation$parallel_workers)
shard_id     <- .parse_arg(cli_args, "--shard-id",     default = 1L)
total_shards <- .parse_arg(cli_args, "--total-shards", default = 1L)
stopifnot(
  "shard-id must be >= 1"        = shard_id >= 1L,
  "shard-id must be <= total-shards" = shard_id <= total_shards,
  "total-shards must be >= 1"    = total_shards >= 1L
)

# Reconfigure future plan if user overrode workers
if (n_workers != CONFIG$simulation$parallel_workers) {
  future::plan(future::multisession, workers = n_workers)
}

cli::cli_h1("Run all sims")
cli::cli_alert_info("n_sims per (scenario, design) = {n_sims}")
cli::cli_alert_info("parallel workers              = {n_workers}")
cli::cli_alert_info("shard                         = {shard_id} / {total_shards}")

# --- load rpact design parameters -------------------------------------------
design_obj <- readRDS(file.path(PATH_SIMS, "design_params.rds"))
efficacy_boundary_final <- design_obj$boundaries_tbl$efficacy_z_boundary[
  nrow(design_obj$boundaries_tbl)
]
interim_event_target <- design_obj$interim_event_target
cli::cli_alert_info("efficacy_boundary_final (z)   = {round(efficacy_boundary_final, 4)}")
cli::cli_alert_info("interim_event_target          = {interim_event_target} events")

# --- build simulation grid ---------------------------------------------------
scenarios_df <- purrr::map_dfr(
  CONFIG$scenarios,
  \(s) tibble::tibble(scenario = s$name, hr_true = s$hr)
)

# Stride-assign sim_ids to shards so each shard gets a near-equal subset.
# sim_id seeds are still derived from the original sim_id (line below), so
# the union of all shard outputs is identical to a single un-sharded run.
all_sim_ids       <- seq_len(n_sims)
this_shard_simids <- all_sim_ids[((all_sim_ids - 1L) %% total_shards) + 1L == shard_id]

grid <- tidyr::expand_grid(
  scenarios_df,
  design = c("fixed", "adaptive"),
  sim_id = this_shard_simids
) |>
  dplyr::mutate(
    seed = CONFIG$simulation$seed + sim_id * 10L + as.integer(factor(design))
  )

cli::cli_alert_info("Total simulations to run      = {nrow(grid)} (this shard)")

# --- dispatcher: route one row to the appropriate simulator ------------------
# Error-row builders (.error_row_fixed / .error_row_adaptive) come from
# R/helpers.R, sourced via R/00_setup.R.
dispatch_one <- function(scenario, hr_true, design, sim_id, seed) {
  tryCatch(
    if (design == "fixed") {
      sim_fixed_trial(
        hr_true                 = hr_true,
        n_per_arm               = CONFIG$trial$max_n / 2L,
        accrual_rate            = CONFIG$trial$accrual_rate_per_month,
        follow_months           = CONFIG$trial$total_followup_months,
        baseline_annual_hazard  = CONFIG$trial$baseline_event_rate,
        alpha_one_sided         = CONFIG$design$alpha_one_sided,
        efficacy_boundary_final = efficacy_boundary_final,
        sim_id                  = sim_id,
        seed                    = seed
      ) |> dplyr::mutate(scenario = scenario, .before = 1)
    } else {
      sim_adaptive_trial(
        hr_true                 = hr_true,
        max_n                   = CONFIG$trial$max_n,
        interim_event_target    = interim_event_target,
        futility_p_hr_lt_0_7    = CONFIG$design$futility_boundary,
        accrual_rate            = CONFIG$trial$accrual_rate_per_month,
        follow_months           = CONFIG$trial$total_followup_months,
        baseline_annual_hazard  = CONFIG$trial$baseline_event_rate,
        alpha_one_sided         = CONFIG$design$alpha_one_sided,
        efficacy_boundary_final = efficacy_boundary_final,
        sim_id                  = sim_id,
        seed                    = seed
      ) |> dplyr::mutate(scenario = scenario, .before = 1)
    },
    error = function(e) {
      msg <- conditionMessage(e)
      cli::cli_alert_danger("sim failed (scenario={scenario}, design={design}, sim_id={sim_id}): {msg}")
      if (design == "fixed") {
        .error_row_fixed(scenario, sim_id, seed, hr_true, msg)
      } else {
        .error_row_adaptive(scenario, sim_id, seed, hr_true, msg)
      }
    }
  )
}

# --- run in parallel ---------------------------------------------------------
tic_step("All simulations")
results <- furrr::future_pmap(
  grid |> dplyr::select(scenario, hr_true, design, sim_id, seed),
  dispatch_one,
  .options  = furrr::furrr_options(seed = TRUE),
  .progress = TRUE
) |> dplyr::bind_rows()
toc_step()

cli::cli_alert_success("Got {nrow(results)} result rows")

# --- surface any per-sim failures -------------------------------------------
n_errors <- sum(!is.na(results$error_msg))
if (n_errors > 0) {
  cli::cli_alert_danger("{n_errors} of {nrow(results)} sims errored (decision == 'error')")
  err_summary <- results |>
    dplyr::filter(!is.na(error_msg)) |>
    dplyr::count(scenario, design, error_msg, name = "n_errors") |>
    dplyr::arrange(dplyr::desc(n_errors))
  print(err_summary, n = Inf)
} else {
  cli::cli_alert_success("0 simulation errors")
}

# --- save one RDS per (scenario, design), suffixed with shard if applicable --
shard_suffix <- if (total_shards > 1L) sprintf("_shard%d", shard_id) else ""
results |>
  dplyr::group_by(scenario, design) |>
  dplyr::group_walk(\(df, key) {
    out_path <- file.path(
      PATH_SIMS,
      sprintf("raw_%s_%s%s.rds", key$scenario, key$design, shard_suffix)
    )
    # group_walk strips the group columns from df; bind them back so the saved
    # tibble is self-describing.
    saveRDS(dplyr::bind_cols(key, df), out_path)
    cli::cli_alert_success("Wrote {nrow(df)} rows -> {.path {out_path}}")
  })

# --- quick summary printout --------------------------------------------------
summary_tbl <- results |>
  dplyr::group_by(scenario, design) |>
  dplyr::summarise(
    n_sims      = dplyr::n(),
    reject_rate = mean(decision == "reject", na.rm = TRUE),
    futility    = mean(decision == "futility_stop", na.rm = TRUE),
    mean_n      = mean(n_used, na.rm = TRUE),
    mean_events = mean(n_events, na.rm = TRUE),
    mean_hr_est = mean(hr_est, na.rm = TRUE),
    .groups     = "drop"
  )
cli::cli_h2("Quick OC preview")
print(summary_tbl, n = Inf)

# --- shut down furrr workers cleanly (avoids end-of-run unserialize noise) ---
# multisession workers can race with the parent R process at script exit,
# producing "Error in unserialize(node$con) : error reading from connection".
# Forcing plan(sequential) here closes the cluster before Rscript terminates.
future::plan(future::sequential)
