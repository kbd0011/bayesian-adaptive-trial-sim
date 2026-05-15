# -----------------------------------------------------------------------------
# Purpose : Sensitivity of operating characteristics to the futility threshold.
#           Sweeps the Bayesian futility rule P(HR < 0.7 | data) < theta across
#           4 candidate thresholds (0.10, 0.15, 0.20, 0.30) and 6 HR scenarios.
#           Smaller theta = more conservative (rarely stop) -> preserves power,
#           costs operational efficiency. Larger theta = more aggressive
#           (frequently stop) -> saves enrollment, may lose power.
# Inputs  : outputs/sims/design_params.rds (interim_event_target, eff bound)
# Outputs : outputs/sims/sensitivity_futility.rds   (raw per-sim results)
#           outputs/tables/futility_sensitivity.csv (aggregated table)
#           outputs/figures/futility_sensitivity.{pdf,png}
# CLI     : --n-sims N  (default 300)
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "03_sim_adaptive.R"))

# --- CLI ---------------------------------------------------------------------
.parse_arg <- function(args, key, default, as = as.integer) {
  i <- which(args == key)
  if (length(i) == 1L && i < length(args)) as(args[i + 1L]) else default
}
n_sims <- .parse_arg(commandArgs(trailingOnly = TRUE), "--n-sims", default = 300L)

# --- thresholds to sweep -----------------------------------------------------
thresholds <- c(0.10, 0.15, 0.20, 0.30)

cli::cli_h1("Futility-threshold sensitivity sweep")
cli::cli_alert_info("n_sims per (scenario, threshold) = {n_sims}")
cli::cli_alert_info("thresholds = {paste(thresholds, collapse = ', ')}")

design_obj   <- readRDS(file.path(PATH_SIMS, "design_params.rds"))
eff_bound    <- design_obj$boundaries_tbl$efficacy_z_boundary[2L]
event_target <- design_obj$interim_event_target

scenarios_df <- purrr::map_dfr(
  CONFIG$scenarios,
  \(s) tibble::tibble(scenario = s$name, hr_true = s$hr)
)

grid <- tidyr::expand_grid(
  scenarios_df,
  futility_threshold = thresholds,
  sim_id = seq_len(n_sims)
) |>
  dplyr::mutate(
    seed = CONFIG$simulation$seed +
           sim_id * 10L +
           as.integer(futility_threshold * 1000L)
  )

cli::cli_alert_info("Total sims to run = {nrow(grid)}")

# --- dispatcher --------------------------------------------------------------
dispatch_sens <- function(scenario, hr_true, futility_threshold, sim_id, seed) {
  tryCatch(
    sim_adaptive_trial(
      hr_true                 = hr_true,
      max_n                   = CONFIG$trial$max_n,
      interim_event_target    = event_target,
      futility_p_hr_lt_0_7    = futility_threshold,
      accrual_rate            = CONFIG$trial$accrual_rate_per_month,
      follow_months           = CONFIG$trial$total_followup_months,
      baseline_annual_hazard  = CONFIG$trial$baseline_event_rate,
      alpha_one_sided         = CONFIG$design$alpha_one_sided,
      efficacy_boundary_final = eff_bound,
      sim_id                  = sim_id,
      seed                    = seed
    ) |>
      dplyr::mutate(scenario = scenario,
                    futility_threshold = futility_threshold,
                    .before = 1),
    error = function(e) {
      tibble::tibble(
        scenario = scenario, futility_threshold = futility_threshold,
        sim_id = sim_id, design = "adaptive", seed = seed, hr_true = hr_true,
        n_used = NA_integer_, n_events = NA_integer_,
        hr_est = NA_real_, log_hr_est = NA_real_, pvalue_logrank = NA_real_,
        decision = "error", runtime_sec = NA_real_,
        interim_decision = NA_character_, interim_n_events = NA_integer_,
        final_alloc_ratio = NA_real_, n_refits = NA_integer_,
        p_hr_lt_07_interim = NA_real_,
        error_msg = conditionMessage(e)
      )
    }
  )
}

# --- run in parallel ---------------------------------------------------------
tic_step("Futility sensitivity sims")
results <- furrr::future_pmap(
  grid |> dplyr::select(scenario, hr_true, futility_threshold, sim_id, seed),
  dispatch_sens,
  .options  = furrr::furrr_options(seed = TRUE),
  .progress = TRUE
) |> dplyr::bind_rows()
toc_step()

cli::cli_alert_success("Got {nrow(results)} sims")
n_errors <- sum(!is.na(results$error_msg))
if (n_errors > 0) {
  cli::cli_alert_danger("{n_errors} errors")
}

# --- aggregate ---------------------------------------------------------------
sens <- results |>
  dplyr::group_by(scenario, hr_true, futility_threshold) |>
  dplyr::summarise(
    n_sims          = dplyr::n(),
    reject_rate     = mean(decision == "reject", na.rm = TRUE),
    p_futility_stop = mean(interim_decision == "stop_futility", na.rm = TRUE),
    e_n             = mean(n_used, na.rm = TRUE),
    e_events        = mean(n_events, na.rm = TRUE),
    mean_hr_est     = mean(hr_est, na.rm = TRUE),
    .groups         = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(hr_true), futility_threshold)

# --- save --------------------------------------------------------------------
rds_path <- file.path(PATH_SIMS, "sensitivity_futility.rds")
csv_path <- file.path(PATH_TBL,  "futility_sensitivity.csv")
saveRDS(list(raw = results, summary = sens), rds_path)
utils::write.csv(sens, csv_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.path {rds_path}}")
cli::cli_alert_success("Wrote {.path {csv_path}}")

cli::cli_h2("Futility-threshold sensitivity ({nrow(sens)} cells)")
print(sens, n = Inf)

# --- plots -------------------------------------------------------------------
# Lock scenario order and label thresholds for readability.
scenario_levels <- sens |>
  dplyr::distinct(scenario, hr_true) |>
  dplyr::arrange(dplyr::desc(hr_true)) |>
  dplyr::pull(scenario)
sens <- sens |>
  dplyr::mutate(
    scenario = factor(scenario, levels = scenario_levels),
    threshold_label = factor(sprintf("θ = %.2f", futility_threshold))
  )

theme_oc <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "top",
    strip.background = ggplot2::element_rect(fill = "grey95", color = NA)
  )

p_fut <- ggplot2::ggplot(
  sens,
  ggplot2::aes(x = futility_threshold, y = p_futility_stop, color = scenario, group = scenario)
) +
  ggplot2::geom_line(alpha = 0.6) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_color_brewer(palette = "Set2") +
  ggplot2::scale_x_continuous(breaks = thresholds) +
  ggplot2::labs(
    title    = "Probability of early futility stop vs threshold",
    subtitle = "θ = upper bound on P(HR < 0.7 | data) for futility decision",
    x        = "Futility threshold θ",
    y        = "Pr(stop for futility)",
    color    = "Scenario"
  ) +
  theme_oc

p_rej <- ggplot2::ggplot(
  sens,
  ggplot2::aes(x = futility_threshold, y = reject_rate, color = scenario, group = scenario)
) +
  ggplot2::geom_hline(yintercept = 0.025, linetype = "dotted", color = "grey50") +
  ggplot2::geom_line(alpha = 0.6) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_color_brewer(palette = "Set2") +
  ggplot2::scale_x_continuous(breaks = thresholds) +
  ggplot2::labs(
    title    = "Rejection rate vs threshold",
    subtitle = "Dotted line = nominal α = 0.025",
    x        = "Futility threshold θ",
    y        = "Pr(reject H0)",
    color    = "Scenario"
  ) +
  theme_oc

p_combined <- patchwork::wrap_plots(p_fut, p_rej, ncol = 1, guides = "collect") &
  ggplot2::theme(legend.position = "top")

pdf_path <- file.path(PATH_FIG, "futility_sensitivity.pdf")
png_path <- file.path(PATH_FIG, "futility_sensitivity.png")
ggplot2::ggsave(pdf_path, p_combined, width = 9, height = 9, dpi = 300)
ggplot2::ggsave(png_path, p_combined, width = 9, height = 9, dpi = 150)
cli::cli_alert_success("Wrote {.path futility_sensitivity.pdf} + .png")

future::plan(future::sequential)
