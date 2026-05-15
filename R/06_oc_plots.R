# -----------------------------------------------------------------------------
# Purpose : Visualize operating characteristics. Produces 4 plots:
#             1. oc_power_curve.pdf   (reject_rate vs scenario, by design)
#             2. oc_expected_n.pdf    (E[N] bar chart by design)
#             3. oc_futility_stop.pdf (P(early stop) for adaptive)
#             4. oc_summary_heatmap.pdf (all metrics tiled, faceted by design)
# Inputs  : outputs/sims/oc_table.rds
# Outputs : 4 PDFs in outputs/figures/ (8x6 in, 300 DPI)
# Author  : Cris Taylor
# Date    : 2026-05-14
# -----------------------------------------------------------------------------

source(here::here("R/00_setup.R"))

tic_step("OC plots")

oc <- readRDS(file.path(PATH_SIMS, "oc_table.rds"))

# Lock scenario order by HR (harmful -> very_strong)
scenario_levels <- oc |>
  dplyr::distinct(scenario, hr_true) |>
  dplyr::arrange(dplyr::desc(hr_true)) |>
  dplyr::pull(scenario)

oc <- oc |>
  dplyr::mutate(
    scenario = factor(scenario, levels = scenario_levels),
    design   = factor(design,   levels = c("fixed", "adaptive"))
  )

theme_oc <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "top",
    strip.background = ggplot2::element_rect(fill = "grey95", color = NA)
  )

# --- 1. Power curve ----------------------------------------------------------
p1 <- ggplot2::ggplot(
  oc, ggplot2::aes(x = scenario, y = reject_rate, color = design, group = design)
) +
  ggplot2::geom_hline(yintercept = 0.025, linetype = "dotted", color = "grey50") +
  ggplot2::geom_hline(yintercept = 0.80,  linetype = "dotted", color = "grey50") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = reject_ci_lo, ymax = reject_ci_hi),
    width = 0.15, position = ggplot2::position_dodge(width = 0.3)
  ) +
  ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0.3)) +
  ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.3), alpha = 0.5) +
  ggplot2::scale_color_brewer(palette = "Set1") +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  ggplot2::annotate("text", x = 0.6, y = 0.06, label = "alpha = 0.025",
                    hjust = 0, size = 3, color = "grey40") +
  ggplot2::annotate("text", x = 0.6, y = 0.82, label = "target power 0.80",
                    hjust = 0, size = 3, color = "grey40") +
  ggplot2::labs(
    title    = "Rejection rate by scenario and design",
    subtitle = "Error bars are exact binomial 95% CIs",
    x        = "Scenario (ordered by true HR, harmful to very strong)",
    y        = "Pr(reject H0)",
    color    = "Design"
  ) +
  theme_oc

ggplot2::ggsave(file.path(PATH_FIG, "oc_power_curve.pdf"), p1,
                width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file.path(PATH_FIG, "oc_power_curve.png"), p1,
                width = 8, height = 6, dpi = 150)
cli::cli_alert_success("Wrote {.path oc_power_curve.pdf} + .png")

# --- 2. Expected N -----------------------------------------------------------
max_n <- CONFIG$trial$max_n
p2 <- ggplot2::ggplot(
  oc, ggplot2::aes(x = scenario, y = e_n, fill = design)
) +
  ggplot2::geom_hline(yintercept = max_n, linetype = "dotted", color = "grey50") +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.0f", e_n)),
    position = ggplot2::position_dodge(width = 0.75),
    vjust = -0.4, size = 3.2
  ) +
  ggplot2::scale_fill_brewer(palette = "Set1") +
  ggplot2::scale_y_continuous(limits = c(0, max_n * 1.1)) +
  ggplot2::annotate("text", x = 0.6, y = max_n + 4,
                    label = sprintf("max N = %d", max_n),
                    hjust = 0, size = 3, color = "grey40") +
  ggplot2::labs(
    title    = "Expected sample size by scenario and design",
    subtitle = "Adaptive design saves N when stopping early for futility",
    x        = "Scenario",
    y        = "E[N enrolled]",
    fill     = "Design"
  ) +
  theme_oc

ggplot2::ggsave(file.path(PATH_FIG, "oc_expected_n.pdf"), p2,
                width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file.path(PATH_FIG, "oc_expected_n.png"), p2,
                width = 8, height = 6, dpi = 150)
cli::cli_alert_success("Wrote {.path oc_expected_n.pdf} + .png")

# --- 3. Futility stop probability (adaptive only) ----------------------------
oc_ad <- oc |> dplyr::filter(design == "adaptive")
p3 <- ggplot2::ggplot(
  oc_ad, ggplot2::aes(x = scenario, y = p_futility_stop, fill = scenario)
) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", p_futility_stop)),
    vjust = -0.4, size = 3.5
  ) +
  ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::scale_y_continuous(limits = c(0, max(oc_ad$p_futility_stop, na.rm = TRUE) * 1.2)) +
  ggplot2::labs(
    title    = "Probability of early stop for futility (adaptive design)",
    subtitle = "Computed at interim using P(HR < 0.7 | data) < 0.20",
    x        = "Scenario",
    y        = "Pr(stop for futility)"
  ) +
  theme_oc

ggplot2::ggsave(file.path(PATH_FIG, "oc_futility_stop.pdf"), p3,
                width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file.path(PATH_FIG, "oc_futility_stop.png"), p3,
                width = 8, height = 6, dpi = 150)
cli::cli_alert_success("Wrote {.path oc_futility_stop.pdf} + .png")

# --- 4. Summary heatmap ------------------------------------------------------
oc_long <- oc |>
  dplyr::select(scenario, design, reject_rate, e_n, e_events, mean_hr_est,
                bias_log_hr, p_futility_stop) |>
  tidyr::pivot_longer(-c(scenario, design), names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    metric = factor(metric, levels = c(
      "reject_rate", "e_n", "e_events",
      "mean_hr_est", "bias_log_hr", "p_futility_stop"
    ), labels = c(
      "Pr(reject)", "E[N]", "E[events]",
      "Mean HR", "Bias (log-HR)", "Pr(futility)"
    )),
    label = ifelse(is.na(value), "", sprintf("%.2f", value))
  )

p4 <- ggplot2::ggplot(
  oc_long, ggplot2::aes(x = scenario, y = metric, fill = value)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2) +
  ggplot2::facet_wrap(~ design, ncol = 1) +
  ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90") +
  ggplot2::labs(
    title    = "OC summary heatmap",
    subtitle = "Cells annotated with values; lighter = higher",
    x        = "Scenario",
    y        = NULL,
    fill     = "Value"
  ) +
  theme_oc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
  )

ggplot2::ggsave(file.path(PATH_FIG, "oc_summary_heatmap.pdf"), p4,
                width = 9, height = 7, dpi = 300)
ggplot2::ggsave(file.path(PATH_FIG, "oc_summary_heatmap.png"), p4,
                width = 9, height = 7, dpi = 150)
cli::cli_alert_success("Wrote {.path oc_summary_heatmap.pdf} + .png")

toc_step()
