# -----------------------------------------------------------------------------
# Shared helpers, sourced by R/00_setup.R so every downstream file has them.
# Currently holds the schema-matching error-row builders for the simulator
# dispatch path (R/04) — exposed here so test-failure-modes.R can verify the
# schema without sourcing the full orchestrator.
# -----------------------------------------------------------------------------

# Failure-path tibble for a fixed-design sim. Matches the column shape of
# sim_fixed_trial() exactly, with NA in the analysis fields and a non-NA
# error_msg so R/05 can detect and surface failures.
.error_row_fixed <- function(scenario, sim_id, seed, hr_true, msg) {
  tibble::tibble(
    scenario       = scenario, sim_id = sim_id, design = "fixed",
    seed           = seed, hr_true = hr_true,
    n_used         = NA_integer_, n_events = NA_integer_,
    hr_est         = NA_real_, log_hr_est = NA_real_,
    pvalue_logrank = NA_real_,
    decision       = "error", runtime_sec = NA_real_,
    error_msg      = msg
  )
}

# Failure-path tibble for an adaptive-design sim. Matches the column shape
# of sim_adaptive_trial() including the adaptive-only fields.
.error_row_adaptive <- function(scenario, sim_id, seed, hr_true, msg) {
  tibble::tibble(
    scenario           = scenario, sim_id = sim_id, design = "adaptive",
    seed               = seed, hr_true = hr_true,
    n_used             = NA_integer_, n_events = NA_integer_,
    hr_est             = NA_real_, log_hr_est = NA_real_,
    pvalue_logrank     = NA_real_,
    decision           = "error", runtime_sec = NA_real_,
    interim_decision   = NA_character_, interim_n_events = NA_integer_,
    final_alloc_ratio  = NA_real_, n_refits = NA_integer_,
    p_hr_lt_07_interim = NA_real_,
    error_msg          = msg
  )
}
