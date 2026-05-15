# -----------------------------------------------------------------------------
# Purpose : Compute group-sequential design parameters (O'Brien-Fleming
#           efficacy + beta-spending futility, one interim at 50% info).
#           Also derive a survival-specific sample size at HR = 0.70 for
#           cross-validation against SAS PROC SEQDESIGN.
# Inputs  : config.yml
# Outputs : outputs/sims/design_params.rds        (rpact design object + extras)
#           outputs/tables/design_boundaries.csv  (tidy stage summary)
#           outputs/figures/alpha_spending.pdf    (alpha-spending function)
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))
suppressPackageStartupMessages(library(rpact))

tic_step("Design parameters")

# --- group-sequential design -------------------------------------------------
design <- rpact::getDesignGroupSequential(
  kMax              = 2L,
  informationRates  = c(CONFIG$design$interim_information_fraction, 1),
  alpha             = CONFIG$design$alpha_one_sided,
  sided             = 1L,
  beta              = 1 - CONFIG$design$target_power,
  typeOfDesign      = "asOF",   # O'Brien-Fleming alpha spending
  typeBetaSpending  = "bsOF"    # O'Brien-Fleming beta spending (futility)
)

cli::cli_h2("rpact group-sequential design")
print(summary(design))

# --- survival sample-size calculation ---------------------------------------
# Convert annual control hazard -> per-month hazard (time unit = months).
lambda2_month <- CONFIG$trial$baseline_event_rate / 12
accrual_time  <- CONFIG$trial$max_n / CONFIG$trial$accrual_rate_per_month
followup_time <- CONFIG$trial$total_followup_months - accrual_time
stopifnot(followup_time > 0)

ss <- rpact::getSampleSizeSurvival(
  design                  = design,
  hazardRatio             = 0.70,
  lambda2                 = lambda2_month,
  accrualTime             = accrual_time,
  followUpTime            = followup_time,
  allocationRatioPlanned  = 1,
  dropoutRate1            = 0,
  dropoutRate2            = 0
)

cli::cli_h2("rpact survival sample size (target HR = 0.70)")
print(summary(ss))

# --- tidy boundaries table ---------------------------------------------------
boundaries_tbl <- tibble::tibble(
  stage                    = seq_len(design$kMax),
  info_fraction            = design$informationRates,
  efficacy_z_boundary      = as.numeric(design$criticalValues),
  futility_z_boundary      = as.numeric(design$futilityBounds[seq_len(design$kMax)]),
  cumulative_alpha_spent   = as.numeric(design$alphaSpent),
  cumulative_beta_spent    = as.numeric(design$betaSpent)
)

readr_path <- file.path(PATH_TBL, "design_boundaries.csv")
utils::write.csv(boundaries_tbl, readr_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.path {readr_path}}")

# --- alpha-spending plot ----------------------------------------------------
fig_path <- file.path(PATH_FIG, "alpha_spending.pdf")
pdf(fig_path, width = 8, height = 6)
print(plot(design, type = 4L))  # type 4 = error-spending function
dev.off()
cli::cli_alert_success("Wrote {.path {fig_path}}")

# --- interim event target (information-time driven) -------------------------
# rpact's `informationRates` are calibrated to the design ALTERNATIVE, so the
# event target is `interim_information_fraction` x expected events under the
# design alternative HR (defaults to 0.70). Closed-form expectation under
# exponential survival with uniform accrual over [0, T_acc] and admin
# censoring at T_tot:
#
#   E[event prob per subject]
#     = 1 - (1/(lambda * T_acc)) * (exp(-lambda*(T_tot-T_acc)) - exp(-lambda*T_tot))
#
# under arm-averaged hazard `lambda_avg = lambda_c * (1 + design_alt_hr) / 2`
# (1:1 allocation pre-interim).
design_alt_hr   <- CONFIG$design$design_alternative_hr
lambda_avg_h1   <- lambda2_month * (1 + design_alt_hr) / 2
T_tot           <- CONFIG$trial$total_followup_months
p_event_h1      <- 1 - (1 / (lambda_avg_h1 * accrual_time)) *
                   (exp(-lambda_avg_h1 * (T_tot - accrual_time)) -
                    exp(-lambda_avg_h1 *  T_tot))
expected_total_events_h1 <- CONFIG$trial$max_n * p_event_h1
interim_event_target     <- ceiling(expected_total_events_h1 *
                                      CONFIG$design$interim_information_fraction)
# Also compute H0 expectation for context in the report.
p_event_h0 <- 1 - (1 / (lambda2_month * accrual_time)) *
              (exp(-lambda2_month * (T_tot - accrual_time)) -
               exp(-lambda2_month *  T_tot))
expected_total_events_h0 <- CONFIG$trial$max_n * p_event_h0

cli::cli_h2("Interim event target ({CONFIG$design$interim_information_fraction*100}% info under H1, HR={design_alt_hr})")
cli::cli_alert_info("E[total events | H1, n={CONFIG$trial$max_n}] = {round(expected_total_events_h1, 1)}")
cli::cli_alert_info("E[total events | H0, n={CONFIG$trial$max_n}] = {round(expected_total_events_h0, 1)}")
cli::cli_alert_info("Interim fires at event count = {interim_event_target}")

# --- persist design object --------------------------------------------------
design_path <- file.path(PATH_SIMS, "design_params.rds")
saveRDS(
  list(
    design                   = design,
    sample_size              = ss,
    boundaries_tbl           = boundaries_tbl,
    lambda2_month            = lambda2_month,
    accrual_time             = accrual_time,
    followup_time            = followup_time,
    target_hr                = 0.70,
    expected_total_events_h0 = expected_total_events_h0,
    expected_total_events_h1 = expected_total_events_h1,
    design_alt_hr            = design_alt_hr,
    interim_event_target     = interim_event_target
  ),
  design_path
)
cli::cli_alert_success("Wrote {.path {design_path}}")

toc_step()
