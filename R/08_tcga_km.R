# -----------------------------------------------------------------------------
# Purpose : Kaplan-Meier analysis of TCGA-BRCA overall survival stratified by
#           hormone-receptor status, with log-rank test and median survival.
# Inputs  : data/tcga_brca_clean.rds
# Outputs : outputs/figures/tcga_km.pdf
#           outputs/tables/tcga_km_summary.csv
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))

tic_step("TCGA KM")

tcga <- readRDS(file.path(PATH_DATA, "tcga_brca_clean.rds"))

# --- Fit KM ------------------------------------------------------------------
fit <- survival::survfit(
  survival::Surv(times, status) ~ hr_status,
  data = tcga
)

# --- ggsurvplot --------------------------------------------------------------
gg <- survminer::ggsurvplot(
  fit,
  data             = tcga,
  risk.table       = TRUE,
  pval             = TRUE,
  pval.method      = TRUE,
  conf.int         = TRUE,
  surv.median.line = "hv",
  palette          = "jco",
  legend.title     = "HR status",
  legend.labs      = levels(factor(tcga$hr_status)),
  xlab             = "Time (days)",
  ylab             = "Overall survival probability",
  ggtheme          = ggplot2::theme_bw(base_size = 11),
  risk.table.height = 0.25
)

# survminer returns a "ggsurvplot" object; print via its built-in method.
fig_path <- file.path(PATH_FIG, "tcga_km.pdf")
pdf(fig_path, width = 8, height = 8)
print(gg)
dev.off()
cli::cli_alert_success("Wrote {.path {fig_path}}")

# --- Summary table: median + 95% CI per stratum + log-rank p ----------------
km_tbl <- summary(fit)$table |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "stratum") |>
  tibble::as_tibble()

lr <- survival::survdiff(survival::Surv(times, status) ~ hr_status, data = tcga)
lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1L)

km_tbl <- km_tbl |>
  dplyr::mutate(
    stratum     = sub("^hr_status=", "", stratum),
    logrank_p   = lr_p,
    median_days = median,
    median_ci   = sprintf("(%g, %g)",
                          `0.95LCL`,
                          `0.95UCL`)
  ) |>
  dplyr::select(stratum, n.start, events, median_days, median_ci, logrank_p)

tbl_path <- file.path(PATH_TBL, "tcga_km_summary.csv")
utils::write.csv(km_tbl, tbl_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.path {tbl_path}}")

cli::cli_h2("TCGA-BRCA KM summary")
print(km_tbl)

toc_step()
