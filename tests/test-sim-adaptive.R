# -----------------------------------------------------------------------------
# Tests for R/03_sim_adaptive.R: sim_adaptive_trial()
# Run via: testthat::test_file("tests/test-sim-adaptive.R")
# Note: This test sources R/03 which compiles a Stan model on first load.
# -----------------------------------------------------------------------------

source(here::here("R", "03_sim_adaptive.R"))

# Boundary needed for the adaptive simulator; load from R/01's output if it
# exists, otherwise use the documented value from the OBF design (z = 1.969
# at final stage for kMax=2, alpha=0.025 one-sided).
design_obj <- tryCatch(
  readRDS(file.path(PATH_SIMS, "design_params.rds")),
  error = function(e) list(
    boundaries_tbl = data.frame(efficacy_z_boundary = c(NA, 1.9686)),
    interim_event_target = 23L
  )
)
eff_bound        <- design_obj$boundaries_tbl$efficacy_z_boundary[2L]
interim_target_t <- design_obj$interim_event_target

test_that("sim_adaptive_trial returns a tibble with the expected columns", {
  res <- sim_adaptive_trial(
    hr_true = 0.70, max_n = 120, interim_event_target = interim_target_t,
    futility_p_hr_lt_0_7 = 0.20, accrual_rate = 12,
    follow_months = 24, baseline_annual_hazard = 0.30,
    alpha_one_sided = 0.025, efficacy_boundary_final = eff_bound,
    sim_id = 1L, seed = 42L
  )
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1)
  expect_equal(res$design, "adaptive")
  expect_true(res$decision %in% c("reject", "accept H0", "futility_stop"))
  expect_true(res$interim_decision %in% c("continue", "stop_futility"))
  expect_true(res$final_alloc_ratio >= 0 && res$final_alloc_ratio <= 1)
})

test_that("sim_adaptive_trial is reproducible given a fixed seed", {
  res_a <- sim_adaptive_trial(
    hr_true = 0.70, max_n = 120, interim_event_target = interim_target_t,
    futility_p_hr_lt_0_7 = 0.20, accrual_rate = 12,
    follow_months = 24, baseline_annual_hazard = 0.30,
    alpha_one_sided = 0.025, efficacy_boundary_final = eff_bound,
    sim_id = 1L, seed = 7777L
  )
  res_b <- sim_adaptive_trial(
    hr_true = 0.70, max_n = 120, interim_event_target = interim_target_t,
    futility_p_hr_lt_0_7 = 0.20, accrual_rate = 12,
    follow_months = 24, baseline_annual_hazard = 0.30,
    alpha_one_sided = 0.025, efficacy_boundary_final = eff_bound,
    sim_id = 1L, seed = 7777L
  )
  expect_equal(res_a$decision, res_b$decision)
  expect_equal(res_a$interim_decision, res_b$interim_decision)
  expect_equal(res_a$n_used, res_b$n_used)
})

test_that(".alloc_prob respects floor/ceiling and applies sqrt damping", {
  expect_equal(.alloc_prob(1.0), 0.8)              # ceiling cap
  expect_equal(.alloc_prob(0.0), 0.2)              # floor cap
  expect_equal(.alloc_prob(0.5), sqrt(0.5))        # mid-range = sqrt
  expect_equal(.alloc_prob(0.25), 0.5)             # sqrt(0.25)
  expect_equal(.alloc_prob(0.04), 0.2)             # sqrt(0.04)=0.2 -> floor
  expect_equal(.alloc_prob(0.64), 0.8)             # sqrt(0.64)=0.8 -> ceiling
})

test_that(".init_adaptive_state returns expected structure", {
  s <- .init_adaptive_state(max_n = 100L)
  expect_length(s$arm, 100L)
  expect_length(s$enrol, 100L)
  expect_length(s$evt_rel, 100L)
  expect_equal(s$alloc_prob_treat, 0.5)
  expect_equal(s$current_time, 0)
  expect_equal(s$n_enrolled, 0L)
  expect_equal(s$n_refits, 0L)
  expect_equal(s$interim_decision, "continue")
  expect_true(is.na(s$n_at_interim))
  expect_false(s$stopped_early)
})

test_that(".observed_events_count counts only events before current_time and follow_months", {
  s <- .init_adaptive_state(5L)
  # 5 subjects all enrolled at time 1, with event times {0.5, 2, 4, 6, 9} relative
  s$enrol[1:5]   <- c(1, 1, 1, 1, 1)
  s$evt_rel[1:5] <- c(0.5, 2, 4, 6, 9)  # calendar event times: 1.5, 3, 5, 7, 10
  s$current_time <- 5
  # Events with cal_time < 5 AND cal_time < follow_months: 1.5 and 3 only.
  expect_equal(.observed_events_count(s, 5L, follow_months = 24), 2L)
  # Tighten follow_months; nothing changes since 24 > 5.
  expect_equal(.observed_events_count(s, 5L, follow_months = 4), 2L)
  # Advance time past follow_months: still bounded by follow_months in calendar.
  s$current_time <- 8
  expect_equal(.observed_events_count(s, 5L, follow_months = 4), 2L)
})

test_that(".snapshot_data correctly censors at current_time", {
  s <- .init_adaptive_state(3L)
  s$arm[1:3]     <- c("control", "treatment", "control")
  s$enrol[1:3]   <- c(0, 1, 2)
  s$evt_rel[1:3] <- c(2, 5, 10)         # cal event times: 2, 6, 12
  s$current_time <- 4                    # max followup: 4, 3, 2
  d <- .snapshot_data(s, 3L)
  expect_equal(d$time,   c(2, 3, 2))     # min(evt_rel, max_followup)
  expect_equal(d$status, c(1L, 0L, 0L))  # only subject 1's event observed (2 < 4)
  expect_equal(d$treat,  c(0L, 1L, 0L))
})

test_that(".enroll_one preserves RNG order vs the original inline loop", {
  # Lockstep: the inline loop calls rexp(accrual), runif(arm), rexp(event).
  # Calling .enroll_one consumes the same 3 RNG draws in the same order.
  s <- .init_adaptive_state(2L)
  set.seed(123L)
  expected <- list(
    t_accrual = rexp(1L, rate = 12),
    u_arm     = runif(1L),
    t_event   = rexp(1L, rate = 0.025)  # control hazard since u_arm here is the assignment
  )
  set.seed(123L)
  s2 <- .enroll_one(s, i = 1L, accrual_rate = 12, lambda_c = 0.025,
                    hr_true = 1.0, follow_months = 24)
  expect_equal(s2$current_time, expected$t_accrual)
  expect_equal(s2$enrol[1L], expected$t_accrual)
  # arm depends on u_arm vs alloc_prob_treat=0.5
  expect_equal(s2$arm[1L],
               if (expected$u_arm < 0.5) "treatment" else "control")
  expect_equal(s2$evt_rel[1L], expected$t_event)
  expect_equal(s2$n_enrolled, 1L)
})

test_that("futility stops under strongly harmful direction (HR=1.5, 10 sims)", {
  # At HR=1.15 the validated futility rate is 48%; HR=1.5 should be higher
  # (~60-70%). On only 10 sims the binomial MC SE is ~0.15, so we require
  # >= 4 stops (a 3-sigma lower bound around the expected rate). A genuinely
  # broken simulator would produce <= 1 stop here.
  res <- purrr::map_dfr(seq_len(10), function(i) {
    sim_adaptive_trial(
      hr_true = 1.50, max_n = 120, interim_event_target = interim_target_t,
      futility_p_hr_lt_0_7 = 0.20, accrual_rate = 12,
      follow_months = 24, baseline_annual_hazard = 0.30,
      alpha_one_sided = 0.025, efficacy_boundary_final = eff_bound,
      sim_id = i, seed = 20000L + i
    )
  })
  expect_gte(sum(res$interim_decision == "stop_futility"), 4)
})
