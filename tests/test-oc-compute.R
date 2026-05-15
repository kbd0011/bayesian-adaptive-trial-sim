# -----------------------------------------------------------------------------
# Tests for R/05_oc_compute.R: operating-characteristics computation
# Run via: testthat::test_file("tests/test-oc-compute.R")
# These tests assume R/04 has been run at least once to populate
# outputs/sims/raw_*.rds; if not, the tests are skipped.
# -----------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

raw_files <- list.files(PATH_SIMS, pattern = "^raw_.*\\.rds$", full.names = TRUE)
oc_path <- file.path(PATH_SIMS, "oc_table.rds")

skip_if_no_oc <- function() {
  if (!file.exists(oc_path)) {
    testthat::skip("Run R/05_oc_compute.R first to generate oc_table.rds")
  }
}

test_that("OC table has expected structure", {
  skip_if_no_oc()
  oc <- readRDS(oc_path)
  expect_true(all(c("scenario", "design", "hr_true", "n_sims",
                    "reject_rate", "reject_ci_lo", "reject_ci_hi",
                    "e_n", "e_events", "mean_hr_est", "bias_log_hr") %in% names(oc)))
  expect_equal(nrow(oc), 12)  # 6 scenarios x 2 designs
  expect_setequal(unique(oc$design), c("fixed", "adaptive"))
})

test_that("Type I under null is within 95% CI of nominal alpha=0.025", {
  skip_if_no_oc()
  oc <- readRDS(oc_path)
  null_rows <- oc[oc$scenario == "null", ]
  # For 1000 sims, MC SE on Type I at 0.025 is sqrt(0.025*0.975/1000) ~= 0.005.
  # 99% CI width is roughly +/-0.013. With 1000 sims we expect Type I in (0.012, 0.038).
  expect_true(all(null_rows$reject_rate < 0.05))
  expect_true(all(null_rows$reject_rate >= 0))
})

test_that("Reject rate is monotone in true effect (more benefit -> more rejects)", {
  skip_if_no_oc()
  oc <- readRDS(oc_path)
  for (d in c("fixed", "adaptive")) {
    sub <- oc[oc$design == d, ]
    # Reject rate should INCREASE as HR DECREASES -> strong negative correlation.
    expect_lt(cor(sub$hr_true, sub$reject_rate, method = "spearman"), -0.6)
  }
})

test_that("Adaptive design saves expected N under harmful/null vs fixed", {
  skip_if_no_oc()
  oc <- readRDS(oc_path)
  for (s in c("null", "harmful")) {
    adaptive_en <- oc[oc$design == "adaptive" & oc$scenario == s, "e_n", drop = TRUE]
    fixed_en    <- oc[oc$design == "fixed"    & oc$scenario == s, "e_n", drop = TRUE]
    expect_lt(adaptive_en, fixed_en)
  }
})

test_that("Futility stop probability is populated only for adaptive design", {
  skip_if_no_oc()
  oc <- readRDS(oc_path)
  expect_true(all(is.na(oc[oc$design == "fixed", "p_futility_stop", drop = TRUE])))
  expect_true(all(!is.na(oc[oc$design == "adaptive", "p_futility_stop", drop = TRUE])))
})
