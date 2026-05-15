# -----------------------------------------------------------------------------
# Purpose : Cox PH regression on TCGA-BRCA overall survival ~ HR_status +
#           age_decade. Proportional hazards diagnostic via Schoenfeld
#           residuals; falls back to stratified Cox as sensitivity if PH
#           assumption is violated for any covariate.
# Inputs  : data/tcga_brca_clean.rds
# Outputs : outputs/tables/tcga_cox_results.csv
#           outputs/figures/tcga_cox_schoenfeld.pdf
#           outputs/sims/tcga_cox_model.rds
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))

tic_step("TCGA Cox PH")

tcga <- readRDS(file.path(PATH_DATA, "tcga_brca_clean.rds"))

# --- Cox PH fit --------------------------------------------------------------
cox <- survival::coxph(
  survival::Surv(times, status) ~ hr_status + age_decade,
  data = tcga,
  ties = "efron"
)
print(summary(cox))

# --- Tidy results ------------------------------------------------------------
tidy_tbl <- broom::tidy(cox, exponentiate = TRUE, conf.int = TRUE) |>
  dplyr::rename(HR = estimate, lower95 = conf.low, upper95 = conf.high)

tbl_path <- file.path(PATH_TBL, "tcga_cox_results.csv")
utils::write.csv(tidy_tbl, tbl_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.path {tbl_path}}")

# --- PH diagnostic: Schoenfeld residuals ------------------------------------
zph <- survival::cox.zph(cox)
cli::cli_h2("Schoenfeld test for PH assumption")
print(zph)

# ggcoxzph returns a list of ggplots (one per covariate + global). Combine.
zph_gg <- survminer::ggcoxzph(zph,
                              ggtheme = ggplot2::theme_bw(base_size = 11))
fig_path <- file.path(PATH_FIG, "tcga_cox_schoenfeld.pdf")
ggplot2::ggsave(fig_path,
                gridExtra::marrangeGrob(grobs = zph_gg, ncol = 1, nrow = length(zph_gg)),
                width = 8, height = 4 * length(zph_gg), dpi = 300,
                limitsize = FALSE)
cli::cli_alert_success("Wrote {.path {fig_path}}")

# --- Sensitivity: stratified Cox if any PH violation -------------------------
violators <- rownames(zph$table)
violator_p <- zph$table[, "p"]
needs_strat <- !is.na(violator_p) &
               violator_p < 0.05 &
               rownames(zph$table) != "GLOBAL"

cox_strat <- NULL
if (any(needs_strat)) {
  # Stratify on the first violating categorical covariate; fall back to age
  # decade strata if age violates.
  first_viol <- rownames(zph$table)[needs_strat][1L]
  cli::cli_alert_warning("PH assumption violated for {first_viol} (p={signif(violator_p[needs_strat][1], 3)}); fitting stratified Cox as sensitivity.")
  strat_form <- as.formula(
    sprintf("Surv(times, status) ~ %s + strata(%s)",
            setdiff(c("hr_status", "age_decade"), first_viol)[1L],
            first_viol)
  )
  cox_strat <- survival::coxph(strat_form, data = tcga, ties = "efron")
  cli::cli_h2("Stratified Cox (sensitivity)")
  print(summary(cox_strat))
}

# --- Persist model object(s) -------------------------------------------------
model_path <- file.path(PATH_SIMS, "tcga_cox_model.rds")
saveRDS(list(cox = cox, cox_strat = cox_strat, zph = zph), model_path)
cli::cli_alert_success("Wrote {.path {model_path}}")

toc_step()
