# -----------------------------------------------------------------------------
# Tests for R/02_sim_fixed.R: sim_fixed_trial()
# Run via: testthat::test_file("tests/test-sim-fixed.R")
# -----------------------------------------------------------------------------

source(here::here("R", "02_sim_fixed.R"))

# Load the rpact boundary; fallback to 1.9686 (OBF kMax=2 alpha=0.025) if not.
eff_bound <- tryCatch(
  readRDS(file.path(PATH_SIMS, "design_params.rds"))$boundaries_tbl$
    efficacy_z_boundary[2L],
  error = function(e) 1.9686
)

test_that("sim_fixed_trial returns a tibble with the expected columns", {
  res <- sim_fixed_trial(
    hr_true = 0.70, n_per_arm = 60,
    accrual_rate = 12, follow_months = 24,
    baseline_annual_hazard = 0.30, alpha_one_sided = 0.025,
    efficacy_boundary_final = eff_bound,
    sim_id = 1L, seed = 42L
  )
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1)
  expect_named(res, c("sim_id", "design", "seed", "hr_true",
                      "n_used", "n_events", "hr_est", "log_hr_est",
                      "pvalue_logrank", "decision", "runtime_sec",
                      "error_msg"))
  expect_equal(res$design, "fixed")
  expect_true(res$decision %in% c("reject", "accept H0"))
})

test_that("sim_fixed_trial is reproducible given a fixed seed", {
  res_a <- sim_fixed_trial(
    hr_true = 0.70, n_per_arm = 60,
    accrual_rate = 12, follow_months = 24,
    baseline_annual_hazard = 0.30, alpha_one_sided = 0.025,
    efficacy_boundary_final = eff_bound,
    sim_id = 1L, seed = 12345L
  )
  res_b <- sim_fixed_trial(
    hr_true = 0.70, n_per_arm = 60,
    accrual_rate = 12, follow_months = 24,
    baseline_annual_hazard = 0.30, alpha_one_sided = 0.025,
    efficacy_boundary_final = eff_bound,
    sim_id = 1L, seed = 12345L
  )
  expect_equal(res_a$decision, res_b$decision)
  expect_equal(res_a$n_used, res_b$n_used)
  expect_equal(res_a$n_events, res_b$n_events)
  expect_equal(res_a$hr_est, res_b$hr_est)
})

test_that("n_used equals 2 * n_per_arm under reasonable params (no truncation)", {
  res <- sim_fixed_trial(
    hr_true = 0.70, n_per_arm = 60,
    accrual_rate = 12, follow_months = 24,
    baseline_annual_hazard = 0.30, alpha_one_sided = 0.025,
    efficacy_boundary_final = eff_bound,
    sim_id = 1L, seed = 99L
  )
  expect_equal(res$n_used, 120)
})

test_that("Type I error under the null is approximately alpha (50 sims)", {
  # True Type I under this design is ~0.021 (validated against 1000-sim run).
  # MC SE on 50 sims = sqrt(0.021 * 0.979 / 50) ~= 0.020; ~3-sigma upper bound
  # is 0.08. A simulator producing Type I >= 0.08 here is meaningfully broken.
  res <- purrr::map_dfr(seq_len(50), function(i) {
    sim_fixed_trial(
      hr_true = 1.00, n_per_arm = 60,
      accrual_rate = 12, follow_months = 24,
      baseline_annual_hazard = 0.30, alpha_one_sided = 0.025,
      efficacy_boundary_final = eff_bound,
      sim_id = i, seed = 1000L + i
    )
  })
  expect_lt(mean(res$decision == "reject"), 0.08)
})
