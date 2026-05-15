# -----------------------------------------------------------------------------
# Integration test for the top-level entry point R/04_run_all_sims.R.
# Actually invokes the orchestrator as a script (not just the underlying
# simulator functions), so YAML-type mismatches, dispatch bugs, schema
# mismatches, etc. are caught here.
#
# Strategy: backup any existing canonical raw_*.rds, run R/04 with a small
# --n-sims, verify file count and schema, then restore the backups.
# -----------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

test_that("R/04_run_all_sims.R end-to-end produces the expected file layout and schema", {
  # Skip if Stan model cache isn't ready — the orchestrator drives R/03 which
  # compiles Stan on first source, adding 30-60 s to the test.
  if (!file.exists(here::here("stan", "exp_survival.rds"))) {
    testthat::skip("Stan cache missing; run R/04 once first to populate stan/exp_survival.rds")
  }
  # Also need design_params.rds for the efficacy boundary + interim target.
  if (!file.exists(file.path(PATH_SIMS, "design_params.rds"))) {
    testthat::skip("design_params.rds missing; run R/01 first")
  }

  # --- backup canonical raw_*.rds so the test doesn't clobber them ---------
  raw_files <- list.files(PATH_SIMS, pattern = "^raw_.*\\.rds$", full.names = TRUE)
  backup_dir <- tempfile("p2_raw_backup_")
  fs::dir_create(backup_dir)
  if (length(raw_files) > 0) {
    fs::file_copy(raw_files, backup_dir, overwrite = TRUE)
    on.exit({
      # Restore from backup; any file the test wrote with the same name will
      # be overwritten with the original 1000-sim version.
      bk <- list.files(backup_dir, pattern = "^raw_.*\\.rds$", full.names = TRUE)
      if (length(bk) > 0) {
        fs::file_copy(bk, PATH_SIMS, overwrite = TRUE)
      }
      fs::dir_delete(backup_dir)
    }, add = TRUE)
  }

  # --- run the orchestrator as a script with a tiny sim count --------------
  rscript <- file.path(R.home("bin"), "Rscript")
  status <- system2(
    rscript,
    args = c(here::here("R", "04_run_all_sims.R"),
             "--n-sims", "5", "--workers", "2"),
    stdout = FALSE, stderr = FALSE
  )
  expect_equal(status, 0,
               info = "R/04 exited non-zero; orchestrator integration broken")

  # --- file layout: 12 raw_*.rds files, one per (scenario, design) ---------
  written <- list.files(PATH_SIMS, pattern = "^raw_.*\\.rds$", full.names = TRUE)
  expect_equal(length(written), 12L,
               info = sprintf("expected 12 raw_*.rds files, got %d", length(written)))

  # --- schema: each file is a tibble with the expected columns -------------
  expected_fixed_cols <- c("scenario", "sim_id", "design", "seed", "hr_true",
                           "n_used", "n_events", "hr_est", "log_hr_est",
                           "pvalue_logrank", "decision", "runtime_sec",
                           "error_msg")
  expected_adaptive_extra <- c("interim_decision", "interim_n_events",
                               "final_alloc_ratio", "n_refits",
                               "p_hr_lt_07_interim")

  for (f in written) {
    d <- readRDS(f)
    expect_true("design" %in% names(d), info = paste("no design col in", f))
    if (unique(d$design) == "fixed") {
      expect_true(all(expected_fixed_cols %in% names(d)),
                  info = paste("fixed-schema columns missing in", basename(f)))
    } else {
      expect_true(all(c(expected_fixed_cols, expected_adaptive_extra) %in% names(d)),
                  info = paste("adaptive-schema columns missing in", basename(f)))
    }
    expect_equal(nrow(d), 5L,
                 info = paste("expected 5 rows in", basename(f), "got", nrow(d)))
  }
})
