# -----------------------------------------------------------------------------
# Purpose : Bayesian Weibull accelerated-failure-time model on TCGA-BRCA.
#           Compares AFT time ratios (transformed to hazard ratios) against
#           Cox PH point estimates as a coherent-model cross-check.
# Inputs  : data/tcga_brca_clean.rds
#           outputs/sims/tcga_cox_model.rds  (for forest-plot overlay)
#           stan/weibull_aft.stan
# Outputs : outputs/sims/tcga_bayes_aft_fit.rds
#           outputs/tables/tcga_bayes_aft_summary.csv
#           outputs/figures/tcga_aft_rhat.pdf
#           outputs/figures/tcga_aft_vs_cox_forest.pdf
#           outputs/figures/tcga_aft_ppc.pdf
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))
suppressPackageStartupMessages({
  library(rstan)
  library(bayesplot)
  library(posterior)
})
options(rstan.auto_write = TRUE)
rstan_options(auto_write = TRUE)

tic_step("TCGA Bayesian Weibull AFT")

tcga <- readRDS(file.path(PATH_DATA, "tcga_brca_clean.rds"))

# --- Design matrix (no intercept; intercept is a Stan parameter) ------------
# Treatment coding: hr_status_HR_plus = 1 if HR+, 0 if HR-.
X <- model.matrix(~ hr_status + age_decade, data = tcga)[, -1L, drop = FALSE]
colnames(X) <- c("hr_status_HR_plus", "age_decade")
K <- ncol(X)
N <- nrow(tcga)
cli::cli_alert_info("Design matrix: {N} rows x {K} cols ({paste(colnames(X), collapse=', ')})")

stan_data <- list(
  N      = N,
  times  = pmax(tcga$times, 0.1),    # Weibull needs strictly positive times
  status = as.integer(tcga$status),
  K      = K,
  X      = X
)

# --- Fit ---------------------------------------------------------------------
stan_model_file <- here::here("stan", "weibull_aft.stan")
stan_mod <- rstan::stan_model(file = stan_model_file, model_name = "weibull_aft")

fit <- rstan::sampling(
  stan_mod, data = stan_data,
  chains = 4, iter = 2000, warmup = 1000,
  cores = CONFIG$simulation$parallel_workers,
  refresh = 200, control = list(adapt_delta = 0.95)
)

# --- Diagnostics: Rhat / ESS ------------------------------------------------
draws <- posterior::as_draws_df(fit)
core_pars <- c("intercept", paste0("beta[", seq_len(K), "]"), "shape",
               paste0("time_ratio[", seq_len(K), "]"))
sum_tbl <- posterior::summarise_draws(
  posterior::subset_draws(draws, variable = core_pars),
  mean, median, sd, ~ quantile(.x, c(0.025, 0.975)), rhat, ess_bulk, ess_tail
)
print(sum_tbl)

cli::cli_alert_info("Max Rhat across core params: {round(max(sum_tbl$rhat, na.rm=TRUE), 4)}")
cli::cli_alert_info("Min bulk ESS: {round(min(sum_tbl$ess_bulk, na.rm=TRUE))}")

# Rhat histogram (use rhat values already computed in sum_tbl)
rhat_vec <- stats::setNames(sum_tbl$rhat, sum_tbl$variable)
rhat_path <- file.path(PATH_FIG, "tcga_aft_rhat.pdf")
ggplot2::ggsave(rhat_path,
  bayesplot::mcmc_rhat(rhat = rhat_vec) +
    ggplot2::labs(title = "Rhat convergence diagnostics",
                  subtitle = "All parameters: Rhat < 1.05 indicates chain convergence") +
    ggplot2::theme_bw(base_size = 11),
  width = 8, height = 5, dpi = 300
)
cli::cli_alert_success("Wrote {.path {rhat_path}}")

# --- Save posterior summary as CSV ------------------------------------------
sum_path <- file.path(PATH_TBL, "tcga_bayes_aft_summary.csv")
utils::write.csv(as.data.frame(sum_tbl), sum_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.path {sum_path}}")

# --- Forest plot: Bayes AFT vs Cox PH on the HR scale -----------------------
# AFT time_ratio: TR>1 means survival extended -> equivalent HR = 1/TR.
beta_draws <- posterior::subset_draws(draws, variable = "beta") |>
  posterior::as_draws_matrix()  # [draws x K]
hr_draws <- exp(-beta_draws)    # HR = 1 / exp(beta) = exp(-beta)
colnames(hr_draws) <- colnames(X)

bayes_summary <- tibble::tibble(
  term     = colnames(hr_draws),
  estimate = apply(hr_draws, 2, median),
  lower95  = apply(hr_draws, 2, quantile, probs = 0.025),
  upper95  = apply(hr_draws, 2, quantile, probs = 0.975),
  source   = "Bayesian Weibull AFT (HR scale)"
)

cox_obj <- readRDS(file.path(PATH_SIMS, "tcga_cox_model.rds"))$cox
cox_summary <- broom::tidy(cox_obj, exponentiate = TRUE, conf.int = TRUE) |>
  dplyr::transmute(
    term      = dplyr::recode(term,
                              "hr_statusHR+" = "hr_status_HR_plus"),
    estimate, lower95 = conf.low, upper95 = conf.high,
    source    = "Cox PH"
  )

forest_df <- dplyr::bind_rows(bayes_summary, cox_summary) |>
  dplyr::mutate(term = factor(term, levels = c("hr_status_HR_plus", "age_decade")))

p_forest <- ggplot2::ggplot(
  forest_df, ggplot2::aes(x = estimate, y = term, color = source, shape = source)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dotted", color = "grey50") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = lower95, xmax = upper95),
    height = 0.15, position = ggplot2::position_dodge(width = 0.5)
  ) +
  ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0.5)) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_color_brewer(palette = "Set1") +
  ggplot2::labs(
    title    = "Hazard ratios: Bayesian Weibull AFT vs Cox PH",
    subtitle = "AFT time ratios inverted (HR = 1/TR) for direct comparison",
    x        = "Hazard ratio (log scale)",
    y        = NULL,
    color    = NULL, shape = NULL
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "top")

forest_path <- file.path(PATH_FIG, "tcga_aft_vs_cox_forest.pdf")
ggplot2::ggsave(forest_path, p_forest, width = 8, height = 5, dpi = 300)
cli::cli_alert_success("Wrote {.path {forest_path}}")

# --- Posterior predictive KM overlay ----------------------------------------
yrep <- rstan::extract(fit, pars = "y_rep", permuted = TRUE)$y_rep  # [draws x N]
# Subsample to keep PP check manageable
keep_draws <- sample(seq_len(nrow(yrep)), min(50, nrow(yrep)))
yrep_sub   <- yrep[keep_draws, , drop = FALSE]

p_ppc <- bayesplot::ppc_km_overlay(
  y        = stan_data$times,
  yrep     = yrep_sub,
  status_y = stan_data$status
) +
  ggplot2::labs(
    title    = "Posterior-predictive KM overlay (50 replicate trials)",
    subtitle = "Bayesian Weibull AFT fit on TCGA-BRCA",
    x        = "Time (days)", y = "Survival probability"
  ) +
  ggplot2::theme_bw(base_size = 11)

ppc_path <- file.path(PATH_FIG, "tcga_aft_ppc.pdf")
ggplot2::ggsave(ppc_path, p_ppc, width = 8, height = 5, dpi = 300)
cli::cli_alert_success("Wrote {.path {ppc_path}}")

# --- Persist fit -------------------------------------------------------------
fit_path <- file.path(PATH_SIMS, "tcga_bayes_aft_fit.rds")
saveRDS(fit, fit_path)
cli::cli_alert_success("Wrote {.path {fit_path}}")

toc_step()
