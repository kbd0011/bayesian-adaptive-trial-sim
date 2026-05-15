# -----------------------------------------------------------------------------
# Purpose : Simulate a single fixed-design Phase II oncology trial with 1:1
#           randomization, exponential survival data, administrative censoring
#           at end of study, and one final analysis (log-rank + Cox PH).
# Inputs  : function arguments (no file I/O).
# Outputs : 1-row tibble: sim_id, design, seed, hr_true, n_used, n_events,
#           hr_est, log_hr_est, pvalue_logrank, decision, runtime_sec.
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))

sim_fixed_trial <- function(hr_true,
                            n_per_arm,
                            accrual_rate,
                            follow_months,
                            baseline_annual_hazard,
                            alpha_one_sided,
                            efficacy_boundary_final,
                            sim_id,
                            seed) {

  stopifnot(
    is.numeric(hr_true), hr_true > 0,
    is.numeric(n_per_arm), n_per_arm >= 5,
    is.numeric(accrual_rate), accrual_rate > 0,
    is.numeric(follow_months), follow_months > 0,
    is.numeric(baseline_annual_hazard), baseline_annual_hazard > 0,
    is.numeric(alpha_one_sided), alpha_one_sided > 0, alpha_one_sided < 0.5,
    is.numeric(efficacy_boundary_final), efficacy_boundary_final > 0
  )

  set.seed(seed)
  t0 <- Sys.time()

  # --- 1:1 randomization (random order to avoid systematic enrollment bias) ---
  n_total <- 2L * as.integer(n_per_arm)
  arm <- sample(rep(c("control", "treatment"), each = n_per_arm))
  arm <- factor(arm, levels = c("control", "treatment"))

  # --- enrollment times: cumsum of Exp(accrual_rate) inter-arrival times -----
  inter_arrival   <- rexp(n_total, rate = accrual_rate)
  enrollment_time <- cumsum(inter_arrival)

  # Drop subjects who would enroll past end of study (rare for max_n << rate*follow)
  enrolled        <- enrollment_time < follow_months
  enrollment_time <- enrollment_time[enrolled]
  arm             <- arm[enrolled]
  n_used          <- sum(enrolled)

  # --- event times: exponential per arm (rates per month) --------------------
  lambda_c <- baseline_annual_hazard / 12
  lambda_t <- lambda_c * hr_true

  event_time <- numeric(n_used)
  is_ctrl    <- arm == "control"
  event_time[ is_ctrl] <- rexp(sum( is_ctrl), rate = lambda_c)
  event_time[!is_ctrl] <- rexp(sum(!is_ctrl), rate = lambda_t)

  # --- administrative censoring at follow_months from study start ------------
  max_followup  <- follow_months - enrollment_time
  observed_time <- pmin(event_time, max_followup)
  status        <- as.integer(event_time < max_followup)
  n_events      <- sum(status)

  # --- analysis: log-rank + Cox PH -------------------------------------------
  surv <- survival::Surv(observed_time, status)

  analysis <- tryCatch({
    sd       <- survival::survdiff(surv ~ arm)
    pv_two   <- 1 - pchisq(sd$chisq, df = 1L)
    cox      <- survival::coxph(surv ~ arm)
    log_hr   <- unname(coef(cox)[1L])
    se       <- unname(sqrt(diag(vcov(cox)))[1L])
    list(pvalue_two_sided = pv_two, log_hr = log_hr, hr = exp(log_hr), se = se)
  }, error = function(e) {
    list(pvalue_two_sided = NA_real_, log_hr = NA_real_, hr = NA_real_, se = NA_real_)
  })

  # Cox z-statistic in "positive = benefit" convention. Reject if z exceeds
  # the rpact OBF final-stage boundary AND HR < 1. Using the same boundary as
  # the adaptive design makes the comparison apples-to-apples.
  z_neg <- if (!is.na(analysis$log_hr) && !is.na(analysis$se) && analysis$se > 0) {
    -analysis$log_hr / analysis$se
  } else NA_real_

  decision <- if (!is.na(z_neg) && z_neg > efficacy_boundary_final &&
                  !is.na(analysis$hr) && analysis$hr < 1) {
    "reject"
  } else {
    "accept H0"
  }

  runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  tibble::tibble(
    sim_id         = sim_id,
    design         = "fixed",
    seed           = seed,
    hr_true        = hr_true,
    n_used         = n_used,
    n_events       = n_events,
    hr_est         = analysis$hr,
    log_hr_est     = analysis$log_hr,
    pvalue_logrank = analysis$pvalue_two_sided,
    decision       = decision,
    runtime_sec    = runtime_sec,
    error_msg      = NA_character_
  )
}

# --- interactive demo: 10 sims under the null + 10 under strong effect -------
if (interactive()) {
  eff_bound <- readRDS(file.path(PATH_SIMS, "design_params.rds"))$boundaries_tbl$
                efficacy_z_boundary[2L]
  cli::cli_h2("Demo: 10 sims under H0 (HR=1.00)")
  demo_null <- purrr::map_dfr(1:10, function(i) {
    sim_fixed_trial(
      hr_true                = 1.00,
      n_per_arm              = CONFIG$trial$max_n / 2,
      accrual_rate           = CONFIG$trial$accrual_rate_per_month,
      follow_months          = CONFIG$trial$total_followup_months,
      baseline_annual_hazard = CONFIG$trial$baseline_event_rate,
      alpha_one_sided        = CONFIG$design$alpha_one_sided,
      efficacy_boundary_final = eff_bound,
      sim_id                 = i,
      seed                   = CONFIG$simulation$seed + i
    )
  })
  print(demo_null)
  cli::cli_alert_info("Null rejection rate: {mean(demo_null$decision == 'reject')}")
}
