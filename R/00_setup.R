# -----------------------------------------------------------------------------
# Purpose : Project-wide setup. Loads libraries, reads config, defines paths,
#           configures rstan + furrr, sets master seed, exposes log/tic helpers.
#           Every other R script begins with: source(here::here("R/00_setup.R"))
# Inputs  : config.yml
# Outputs : CONFIG list in calling env; PATH_* path constants; logging helpers.
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(fs)
  library(cli)
  library(tictoc)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(furrr)
  library(survival)
  library(survminer)
})
# Heavy Bayesian libs (rstan, bayesplot, posterior) are loaded only by the
# files that need them (R/03, R/10, R/11). This shaves ~15-20 s of startup
# off everything else (tests, R/01, R/05, R/06, R/07).

# --- config ------------------------------------------------------------------
CONFIG <- yaml::read_yaml(here::here("config.yml"))

# yaml::read_yaml usually parses types correctly in R, but be explicit so
# downstream arithmetic / RNG seeding never sees a stray character or list.
# Scientific notation (e.g., 2e-4) and any string-quoted number is forced
# to double here.
CONFIG$trial$accrual_rate_per_month <- as.numeric(CONFIG$trial$accrual_rate_per_month)
CONFIG$trial$total_followup_months  <- as.numeric(CONFIG$trial$total_followup_months)
CONFIG$trial$max_n                  <- as.integer(CONFIG$trial$max_n)
CONFIG$trial$baseline_event_rate    <- as.numeric(CONFIG$trial$baseline_event_rate)

CONFIG$design$alpha_one_sided              <- as.numeric(CONFIG$design$alpha_one_sided)
CONFIG$design$target_power                 <- as.numeric(CONFIG$design$target_power)
CONFIG$design$interim_information_fraction <- as.numeric(CONFIG$design$interim_information_fraction)
CONFIG$design$design_alternative_hr        <- as.numeric(CONFIG$design$design_alternative_hr)
CONFIG$design$futility_boundary            <- as.numeric(CONFIG$design$futility_boundary)

CONFIG$simulation$n_sims_full       <- as.integer(CONFIG$simulation$n_sims_full)
CONFIG$simulation$n_sims_ci         <- as.integer(CONFIG$simulation$n_sims_ci)
CONFIG$simulation$seed              <- as.integer(CONFIG$simulation$seed)
CONFIG$simulation$parallel_workers  <- as.integer(CONFIG$simulation$parallel_workers)

# Sanity bounds — fail fast if config is nonsensical.
stopifnot(
  CONFIG$design$alpha_one_sided > 0,
  CONFIG$design$alpha_one_sided < 0.5,
  CONFIG$design$target_power > 0,
  CONFIG$design$target_power < 1,
  CONFIG$design$interim_information_fraction > 0,
  CONFIG$design$interim_information_fraction < 1,
  CONFIG$design$futility_boundary > 0,
  CONFIG$design$futility_boundary < 1,
  CONFIG$design$design_alternative_hr > 0,
  CONFIG$trial$max_n >= 20L,
  CONFIG$simulation$parallel_workers >= 1L
)

# --- paths -------------------------------------------------------------------
PATH_SIMS       <- here::here("outputs", "sims")
PATH_FIG        <- here::here("outputs", "figures")
PATH_TBL        <- here::here("outputs", "tables")
PATH_DATA       <- here::here("data")
PATH_STAN_CACHE <- here::here("outputs", "sims", ".stan_cache")

fs::dir_create(c(PATH_SIMS, PATH_FIG, PATH_TBL, PATH_DATA, PATH_STAN_CACHE))

# --- parallel ----------------------------------------------------------------
options(mc.cores = CONFIG$simulation$parallel_workers)
# rstan auto-write cache is enabled in R/03 and R/10 (the only callers).
future::plan(future::multisession, workers = CONFIG$simulation$parallel_workers)

# --- master seed -------------------------------------------------------------
set.seed(CONFIG$simulation$seed)

# --- helpers -----------------------------------------------------------------
log_step <- function(msg) {
  cli::cli_alert_info("{format(Sys.time(), '%H:%M:%S')} {msg}")
}

tic_step <- function(msg) {
  tictoc::tic(msg)
  log_step(paste("START:", msg))
}

toc_step <- function() {
  res <- tictoc::toc(quiet = TRUE)
  elapsed <- round(res$toc - res$tic, 2)
  log_step(paste("DONE :", res$msg, "—", elapsed, "s"))
  invisible(elapsed)
}

# --- shared helpers ----------------------------------------------------------
source(here::here("R", "helpers.R"))

# Setup-complete announcement is opt-in (interactive sessions or
# `Sys.setenv(P2_VERBOSE_SETUP = "1")` for debugging). Suppressing this in
# CI/test runs removes ~15 redundant lines of log noise per pipeline run.
if (interactive() || nzchar(Sys.getenv("P2_VERBOSE_SETUP"))) {
  cli::cli_alert_success("Setup complete.")
}
