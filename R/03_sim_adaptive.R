# -----------------------------------------------------------------------------
# Purpose : Simulate a single adaptive Phase II oncology trial with
#             - 1:1 randomization until interim,
#             - Bayesian exponential survival fit at interim,
#             - Futility stop if P(HR < 0.7 | data) < threshold,
#             - Response-adaptive randomization (RAR) post-interim, refit
#               every 20 enrollees,
#             - Final analysis: Cox PH compared to rpact OBF z-boundary.
# Inputs  : function arguments + cached/compiled Stan model.
# Outputs : 1-row tibble; columns mirror sim_fixed_trial() plus
#           interim_decision, interim_n_events, final_alloc_ratio,
#           n_refits, p_hr_lt_07_interim.
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))
suppressPackageStartupMessages(library(rstan))
options(rstan.auto_write = TRUE)
rstan_options(auto_write = TRUE)

# --- Compile (or load cached) Stan model at module load ----------------------
# rstan_options(auto_write = TRUE) [set in 00_setup.R] writes the compiled
# binary next to the .stan source; subsequent sourcings pick up the cached
# .rds via stan_model()'s built-in hash-based cache check.
EXP_SURV_STAN_PATH <- here::here("stan", "exp_survival.stan")
stan_exp_model <- rstan::stan_model(file = EXP_SURV_STAN_PATH,
                                    model_name = "exp_survival")

# --- Fit + extract posterior (one batch) -------------------------------------
.fit_exp_post <- function(obs_time, status, treat,
                          n_chains = 2L, n_iter = 1000L, n_warmup = 500L) {
  # Stan rejects exact zeros for exponential; clamp at a tiny positive number
  obs_time_pos <- pmax(obs_time, 1e-4)
  stan_data <- list(
    N     = length(obs_time_pos),
    time  = obs_time_pos,
    evt   = as.integer(status),
    treat = as.integer(treat)
  )
  fit <- rstan::sampling(
    stan_exp_model, data = stan_data,
    chains = n_chains, iter = n_iter, warmup = n_warmup,
    refresh = 0, show_messages = FALSE, verbose = FALSE,
    control = list(adapt_delta = 0.9)
  )
  log_hr_post <- rstan::extract(fit, "log_hr", permuted = TRUE)$log_hr
  hr_post     <- exp(log_hr_post)
  list(
    log_hr_mean    = mean(log_hr_post),
    hr_median      = stats::median(hr_post),
    p_hr_lt_07     = mean(hr_post < 0.7),
    p_treat_better = mean(hr_post < 1.0)
  )
}

# --- Allocation update rule (Thompson-style with damping + caps) ------------
.alloc_prob <- function(p_treat_better, floor = 0.2, ceiling = 0.8) {
  pmin(ceiling, pmax(floor, sqrt(p_treat_better)))
}

# --- State container for the per-trial simulation loop -----------------------
# Holds all mutable per-trial state in a single named list so that the
# enrollment / interim / refit steps can be extracted into self-contained
# helpers below. Keeping these as helpers (vs an inline loop) lets us
# unit-test the pure parts in tests/test-sim-adaptive.R.
.init_adaptive_state <- function(max_n) {
  list(
    arm                = character(max_n),
    enrol              = numeric(max_n),
    evt_rel            = numeric(max_n),         # event time relative to enrollment
    alloc_prob_treat   = 0.5,                    # initial 1:1 RAR
    current_time       = 0,
    n_enrolled         = 0L,
    n_refits           = 0L,
    interim_decision   = "continue",
    interim_n_events   = NA_integer_,
    p_hr_lt_07_interim = NA_real_,
    interim_hr_median  = NA_real_,
    interim_log_hr     = NA_real_,
    stopped_early      = FALSE,
    n_at_interim       = NA_integer_,
    next_refit_at      = NA_integer_
  )
}

# Enroll one subject. Calls rexp(accrual), then (if accrual is still open)
# runif(arm-assignment), then rexp(event-time). RNG order matches the
# original inline loop bit-for-bit. If accrual is past follow-up, returns
# state unchanged except for the advanced current_time (so caller can break).
.enroll_one <- function(state, i, accrual_rate, lambda_c, hr_true, follow_months) {
  state$current_time <- state$current_time + rexp(1L, rate = accrual_rate)
  if (state$current_time >= follow_months) return(state)

  state$enrol[i] <- state$current_time
  state$arm[i]   <- if (runif(1L) < state$alloc_prob_treat) "treatment" else "control"
  lam <- if (state$arm[i] == "treatment") lambda_c * hr_true else lambda_c
  state$evt_rel[i] <- rexp(1L, rate = lam)
  state$n_enrolled <- i
  state
}

# Count observed events at the current calendar time (an event is observed
# iff its calendar time precedes both `current_time` and `follow_months`).
# Pure function; trivially unit-testable.
.observed_events_count <- function(state, i, follow_months) {
  idx <- seq_len(i)
  cal_evt <- state$enrol[idx] + state$evt_rel[idx]
  sum(cal_evt < state$current_time & cal_evt < follow_months)
}

# Snapshot of observed data at the current calendar time, suitable for
# the Bayesian fit. Pure function.
.snapshot_data <- function(state, i) {
  idx <- seq_len(i)
  max_followup <- state$current_time - state$enrol[idx]
  list(
    time   = pmin(state$evt_rel[idx], max_followup),
    status = as.integer(state$evt_rel[idx] < max_followup),
    treat  = as.integer(state$arm[idx] == "treatment")
  )
}

# Fire the interim analysis: Bayesian fit + futility decision + (if continue)
# RAR update + schedule next refit. Returns updated state with stopped_early
# set if the futility rule triggered.
.do_interim <- function(state, i, futility_p_hr_lt_0_7, rar_refit_interval) {
  d <- .snapshot_data(state, i)
  post <- tryCatch(
    .fit_exp_post(d$time, d$status, d$treat),
    error = function(e) {
      warning("Stan fit failed at interim (i=", i, "): ", conditionMessage(e))
      list(log_hr_mean = NA, hr_median = NA, p_hr_lt_07 = NA, p_treat_better = 0.5)
    }
  )
  state$n_refits           <- state$n_refits + 1L
  state$n_at_interim       <- i
  state$interim_n_events   <- sum(d$status)
  state$p_hr_lt_07_interim <- post$p_hr_lt_07
  state$interim_hr_median  <- post$hr_median
  state$interim_log_hr     <- post$log_hr_mean

  if (!is.na(post$p_hr_lt_07) &&
      post$p_hr_lt_07 < futility_p_hr_lt_0_7) {
    state$interim_decision <- "stop_futility"
    state$stopped_early    <- TRUE
  } else {
    state$alloc_prob_treat <- .alloc_prob(post$p_treat_better)
    state$next_refit_at    <- i + rar_refit_interval
  }
  state
}

# Post-interim RAR refit: Bayesian fit on all accrued data, update allocation,
# schedule next refit. Returns updated state.
.do_rar_refit <- function(state, i, rar_refit_interval) {
  d <- .snapshot_data(state, i)
  post <- tryCatch(
    .fit_exp_post(d$time, d$status, d$treat),
    error = function(e) {
      warning("Stan fit failed at refit (i=", i, "): ", conditionMessage(e))
      list(p_treat_better = state$alloc_prob_treat)
    }
  )
  state$n_refits         <- state$n_refits + 1L
  state$alloc_prob_treat <- .alloc_prob(post$p_treat_better)
  state$next_refit_at    <- i + rar_refit_interval
  state
}

# --- Main simulator ----------------------------------------------------------
sim_adaptive_trial <- function(hr_true,
                               max_n,
                               interim_event_target,
                               futility_p_hr_lt_0_7,
                               accrual_rate,
                               follow_months,
                               baseline_annual_hazard,
                               alpha_one_sided,
                               efficacy_boundary_final,
                               sim_id,
                               seed,
                               rar_refit_interval = 20L) {

  stopifnot(
    is.numeric(hr_true), hr_true > 0,
    is.numeric(max_n), max_n >= 20,
    is.numeric(interim_event_target), interim_event_target >= 1,
    is.numeric(futility_p_hr_lt_0_7), futility_p_hr_lt_0_7 > 0, futility_p_hr_lt_0_7 < 1,
    is.numeric(efficacy_boundary_final), efficacy_boundary_final > 0
  )

  set.seed(seed)
  t0 <- Sys.time()

  lambda_c <- baseline_annual_hazard / 12
  state    <- .init_adaptive_state(max_n)

  # --- Enrollment loop (orchestrator) ----------------------------------------
  for (i in seq_len(max_n)) {
    state <- .enroll_one(state, i, accrual_rate, lambda_c, hr_true, follow_months)
    if (state$current_time >= follow_months) break

    if (is.na(state$n_at_interim)) {
      if (.observed_events_count(state, i, follow_months) >= interim_event_target) {
        state <- .do_interim(state, i, futility_p_hr_lt_0_7, rar_refit_interval)
        if (state$stopped_early) break
      }
    } else if (!is.na(state$next_refit_at) &&
               i == state$next_refit_at && i < max_n) {
      state <- .do_rar_refit(state, i, rar_refit_interval)
    }
  }

  # Unpack state for the final-analysis section (kept as locals to minimise
  # downstream diff vs the pre-refactor code path).
  arm_vec            <- state$arm
  enrol_vec          <- state$enrol
  evt_time_rel_vec   <- state$evt_rel
  n_enrolled         <- state$n_enrolled
  n_refits           <- state$n_refits
  interim_decision   <- state$interim_decision
  interim_n_events   <- state$interim_n_events
  p_hr_lt_07_interim <- state$p_hr_lt_07_interim
  interim_hr_median  <- state$interim_hr_median
  interim_log_hr     <- state$interim_log_hr
  stopped_early      <- state$stopped_early
  current_time       <- state$current_time

  # --- Final analysis: assemble dataset at end of follow-up ------------------
  idx_final <- seq_len(n_enrolled)
  enrol_f   <- enrol_vec[idx_final]
  evt_rel_f <- evt_time_rel_vec[idx_final]
  arm_f     <- factor(arm_vec[idx_final], levels = c("control", "treatment"))

  final_alloc_ratio <- mean(arm_f == "treatment")

  if (stopped_early) {
    # Cox on tiny interim data is unstable (separation -> huge |log_hr|). Use
    # the Bayesian posterior that drove the stop decision as the reported HR.
    end_time           <- current_time
    max_followup_final <- end_time - enrol_f
    status_final       <- as.integer(evt_rel_f < max_followup_final)
    n_events           <- sum(status_final)
    cox_res <- list(
      log_hr = interim_log_hr,
      hr     = interim_hr_median,
      se     = NA_real_,
      pvalue = NA_real_
    )
  } else {
    # Trial completed: full Cox PH at end of follow-up
    end_time           <- follow_months
    max_followup_final <- end_time - enrol_f
    observed_time      <- pmin(evt_rel_f, pmax(max_followup_final, 0))
    status_final       <- as.integer(evt_rel_f < max_followup_final)
    n_events           <- sum(status_final)

    cox_res <- tryCatch({
      surv <- survival::Surv(observed_time, status_final)
      cox  <- survival::coxph(surv ~ arm_f)
      log_hr <- unname(coef(cox)[1L])
      se     <- sqrt(diag(vcov(cox)))[1L]
      sd     <- survival::survdiff(surv ~ arm_f)
      pv_two <- 1 - pchisq(sd$chisq, df = 1L)
      list(log_hr = log_hr, hr = exp(log_hr), se = unname(se), pvalue = pv_two)
    }, error = function(e) {
      list(log_hr = NA_real_, hr = NA_real_, se = NA_real_, pvalue = NA_real_)
    })
  }

  # z_neg: positive = treatment benefit. Reject H0 if z_neg > boundary.
  z_neg <- if (!is.na(cox_res$log_hr) && !is.na(cox_res$se) && cox_res$se > 0) {
    -cox_res$log_hr / cox_res$se
  } else NA_real_

  decision <- if (stopped_early) {
    "futility_stop"
  } else if (!is.na(z_neg) && z_neg > efficacy_boundary_final && cox_res$hr < 1) {
    "reject"
  } else {
    "accept H0"
  }

  runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  tibble::tibble(
    sim_id              = sim_id,
    design              = "adaptive",
    seed                = seed,
    hr_true             = hr_true,
    n_used              = n_enrolled,
    n_events            = n_events,
    hr_est              = cox_res$hr,
    log_hr_est          = cox_res$log_hr,
    pvalue_logrank      = cox_res$pvalue,
    decision            = decision,
    runtime_sec         = runtime_sec,
    interim_decision    = interim_decision,
    interim_n_events    = interim_n_events,
    final_alloc_ratio   = final_alloc_ratio,
    n_refits            = n_refits,
    p_hr_lt_07_interim  = p_hr_lt_07_interim,
    error_msg           = NA_character_
  )
}

# --- interactive demo --------------------------------------------------------
if (interactive()) {
  design_obj  <- readRDS(file.path(PATH_SIMS, "design_params.rds"))
  eff_bound   <- design_obj$boundaries_tbl$efficacy_z_boundary[2L]
  event_target <- design_obj$interim_event_target

  cli::cli_h2("Demo: 5 adaptive sims under HR=1.00 (interim at {event_target} events)")
  demo_null <- purrr::map_dfr(1:5, function(i) {
    sim_adaptive_trial(
      hr_true                = 1.00,
      max_n                  = CONFIG$trial$max_n,
      interim_event_target   = event_target,
      futility_p_hr_lt_0_7   = CONFIG$design$futility_boundary,
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
}
