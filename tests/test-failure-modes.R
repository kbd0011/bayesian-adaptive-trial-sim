# -----------------------------------------------------------------------------
# Tests for the failure-mode paths in the simulators and dispatcher.
# Covers (a) input validation (stopifnot triggers), (b) error-row builders
# return a schema matching the success path.
# -----------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))   # also sources R/helpers.R
source(here::here("R", "02_sim_fixed.R"))

test_that("sim_fixed_trial input validation rejects nonsense args", {
  expect_error(
    sim_fixed_trial(
      hr_true = -0.5, n_per_arm = 60, accrual_rate = 12,
      follow_months = 24, baseline_annual_hazard = 0.30,
      alpha_one_sided = 0.025, efficacy_boundary_final = 1.97,
      sim_id = 1L, seed = 1L
    ),
    regexp = "hr_true"
  )
  expect_error(
    sim_fixed_trial(
      hr_true = 0.7, n_per_arm = 2, accrual_rate = 12,
      follow_months = 24, baseline_annual_hazard = 0.30,
      alpha_one_sided = 0.025, efficacy_boundary_final = 1.97,
      sim_id = 1L, seed = 1L
    ),
    regexp = "n_per_arm"
  )
  expect_error(
    sim_fixed_trial(
      hr_true = 0.7, n_per_arm = 60, accrual_rate = 12,
      follow_months = 24, baseline_annual_hazard = 0.30,
      alpha_one_sided = 0.025, efficacy_boundary_final = -1,
      sim_id = 1L, seed = 1L
    ),
    regexp = "efficacy_boundary_final"
  )
})

test_that(".error_row_fixed schema matches sim_fixed_trial success schema (sans scenario)", {
  # Success path schema (with scenario stripped, since dispatch_one adds it)
  ok <- sim_fixed_trial(
    hr_true = 0.7, n_per_arm = 60, accrual_rate = 12, follow_months = 24,
    baseline_annual_hazard = 0.30, alpha_one_sided = 0.025,
    efficacy_boundary_final = 1.97, sim_id = 1L, seed = 1L
  )
  err <- .error_row_fixed(scenario = "x", sim_id = 1L, seed = 1L,
                          hr_true = 0.7, msg = "boom")
  # Error row has scenario at position 1, then the same downstream columns
  # as the success path.
  expect_setequal(names(err)[-1L], names(ok))
  expect_equal(err$decision, "error")
  expect_equal(err$error_msg, "boom")
  expect_true(all(is.na(err$hr_est)))
  expect_true(all(is.na(err$n_used)))
})

test_that(".error_row_adaptive schema matches sim_adaptive_trial success schema (sans scenario)", {
  # Build the adaptive success schema without actually running Stan: just
  # check that the columns line up via a minimal source.
  source(here::here("R", "03_sim_adaptive.R"))
  design <- tryCatch(readRDS(file.path(PATH_SIMS, "design_params.rds")),
                     error = function(e) NULL)
  skip_if(is.null(design), "design_params.rds missing (run R/01 first)")
  ok <- sim_adaptive_trial(
    hr_true = 0.7, max_n = 120,
    interim_event_target = design$interim_event_target,
    futility_p_hr_lt_0_7 = 0.20, accrual_rate = 12,
    follow_months = 24, baseline_annual_hazard = 0.30,
    alpha_one_sided = 0.025,
    efficacy_boundary_final = design$boundaries_tbl$efficacy_z_boundary[2L],
    sim_id = 1L, seed = 1L
  )
  err <- .error_row_adaptive(scenario = "x", sim_id = 1L, seed = 1L,
                             hr_true = 0.7, msg = "boom")
  expect_setequal(names(err)[-1L], names(ok))
  expect_equal(err$decision, "error")
  expect_equal(err$error_msg, "boom")
})

test_that("error rows have the same column count as success rows + scenario + error_msg", {
  ok <- sim_fixed_trial(
    hr_true = 0.7, n_per_arm = 60, accrual_rate = 12, follow_months = 24,
    baseline_annual_hazard = 0.30, alpha_one_sided = 0.025,
    efficacy_boundary_final = 1.97, sim_id = 1L, seed = 1L
  )
  err <- .error_row_fixed("x", 1L, 1L, 0.7, "boom")
  # success has N columns (incl. error_msg=NA); error has N+1 (incl. scenario).
  expect_equal(ncol(err), ncol(ok) + 1L)
  # bind_rows should produce a single shape with no NA-padded extra columns.
  combined <- dplyr::bind_rows(err, dplyr::mutate(ok, scenario = "y", .before = 1))
  expect_equal(ncol(combined), ncol(err))
})
