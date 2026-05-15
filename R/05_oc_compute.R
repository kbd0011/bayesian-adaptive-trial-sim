# -----------------------------------------------------------------------------
# Purpose : Aggregate raw per-trial simulation results into operating
#           characteristics per (scenario, design): Type I / power with 95% CI,
#           expected N, expected events, P(futility stop), mean HR, bias on
#           log-HR scale.
# Inputs  : outputs/sims/raw_<scenario>_<design>.rds (one per cell, 12 total)
# Outputs : outputs/sims/oc_table.rds  (tibble for downstream code)
#           outputs/tables/oc_table.csv (human-readable)
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))

tic_step("OC compute")

# --- load raw sims -----------------------------------------------------------
raw_files <- list.files(PATH_SIMS, pattern = "^raw_.*\\.rds$", full.names = TRUE)
stopifnot("No raw_*.rds files found; run R/04 first" = length(raw_files) > 0)

# Parse scenario + design from filename (design is the LAST segment, allowing
# scenarios with underscores like "very_strong_effect"). Strips an optional
# trailing "_shard<N>" suffix produced by matrix-sharded R/04 runs.
.parse_filename <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))           # "raw_..._<design>[_shard<N>]"
  stem <- sub("^raw_", "", stem)                              # "..._<design>[_shard<N>]"
  stem <- sub("_shard[0-9]+$", "", stem)                      # "..._<design>"
  parts <- strsplit(stem, "_")[[1]]
  list(
    design   = parts[length(parts)],
    scenario = paste(parts[-length(parts)], collapse = "_")
  )
}

raw <- purrr::map_dfr(raw_files, function(f) {
  d  <- readRDS(f)
  fp <- .parse_filename(f)
  # Only inject the columns if missing (be tolerant to both old and fixed R/04 outputs)
  if (!"scenario" %in% names(d)) d$scenario <- fp$scenario
  if (!"design"   %in% names(d)) d$design   <- fp$design
  d
})
cli::cli_alert_info("Loaded {nrow(raw)} sims from {length(raw_files)} files")

# Surface any errored sims so they aren't silently smoothed by na.rm later.
if ("error_msg" %in% names(raw)) {
  n_errors <- sum(!is.na(raw$error_msg))
  if (n_errors > 0) {
    cli::cli_alert_danger(
      "{n_errors} of {nrow(raw)} loaded sims have non-NA error_msg; their analysis fields are NA and will not contribute to OC averages."
    )
  }
}

# Guarantee columns adaptive-only fields exist (NA-filled for fixed)
adaptive_cols <- c("interim_decision", "interim_n_events",
                   "final_alloc_ratio", "n_refits", "p_hr_lt_07_interim")
for (col in adaptive_cols) {
  if (!col %in% names(raw)) raw[[col]] <- NA
}

# --- binom CI helper (vectorized over (x, n)) --------------------------------
.binom_ci <- function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  ci <- stats::binom.test(x, n)$conf.int
  c(lo = ci[1], hi = ci[2])
}

# --- compute OC per (scenario, design) ---------------------------------------
# MCSE columns quantify Monte Carlo precision: binomial proportions use the
# closed-form SE sqrt(p(1-p)/n); continuous summaries use SD/sqrt(n).
oc <- raw |>
  dplyr::group_by(scenario, design, hr_true) |>
  dplyr::summarise(
    n_sims          = dplyr::n(),
    n_reject        = sum(decision == "reject", na.rm = TRUE),
    reject_rate     = n_reject / n_sims,
    mcse_reject     = sqrt(reject_rate * (1 - reject_rate) / n_sims),
    reject_ci_lo    = .binom_ci(n_reject, n_sims)["lo"],
    reject_ci_hi    = .binom_ci(n_reject, n_sims)["hi"],
    p_futility_stop = if (design[1L] == "adaptive") {
                        mean(interim_decision == "stop_futility", na.rm = TRUE)
                      } else NA_real_,
    mcse_p_futility = if (design[1L] == "adaptive") {
                        sqrt(p_futility_stop * (1 - p_futility_stop) / n_sims)
                      } else NA_real_,
    e_n             = mean(n_used,   na.rm = TRUE),
    mcse_e_n        = stats::sd(n_used,   na.rm = TRUE) / sqrt(n_sims),
    e_events        = mean(n_events, na.rm = TRUE),
    mcse_e_events   = stats::sd(n_events, na.rm = TRUE) / sqrt(n_sims),
    mean_hr_est     = mean(hr_est,   na.rm = TRUE),
    median_hr_est   = stats::median(hr_est, na.rm = TRUE),
    bias_log_hr     = mean(log(hr_est[hr_est > 0]), na.rm = TRUE) - log(hr_true[1L]),
    mean_alloc_ratio = mean(final_alloc_ratio, na.rm = TRUE),
    .groups         = "drop"
  ) |>
  dplyr::arrange(design, dplyr::desc(hr_true))

# Tag the alpha vs power column for the CSV
oc <- oc |>
  dplyr::mutate(
    metric_label = ifelse(scenario == "null", "Type I error", "Power"),
    .after = design
  )

# --- save --------------------------------------------------------------------
rds_path <- file.path(PATH_SIMS, "oc_table.rds")
csv_path <- file.path(PATH_TBL,  "oc_table.csv")

saveRDS(oc, rds_path)
utils::write.csv(oc, csv_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.path {rds_path}}")
cli::cli_alert_success("Wrote {.path {csv_path}}")

# --- console preview ---------------------------------------------------------
cli::cli_h2("Operating Characteristics ({nrow(oc)} cells)")
oc_print <- oc |>
  dplyr::transmute(
    scenario, design,
    hr_true,
    rate    = sprintf("%.3f ± %.4f", reject_rate, mcse_reject),
    rate_ci = sprintf("[%.3f, %.3f]", reject_ci_lo, reject_ci_hi),
    fut     = ifelse(is.na(p_futility_stop), "-", sprintf("%.2f", p_futility_stop)),
    eN      = sprintf("%.1f ± %.2f", e_n, mcse_e_n),
    eE      = sprintf("%.1f", e_events),
    HR      = sprintf("%.3f", mean_hr_est),
    bias    = sprintf("%+.3f", bias_log_hr)
  )
print(oc_print, n = Inf)
cli::cli_alert_info("Reported as mean ± MCSE; CIs are exact binomial 95%")

toc_step()
