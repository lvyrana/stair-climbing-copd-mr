suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(svglite)
  library(ragg)
})

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))] |> sub("^--file=", "", x = _)), ".."))
source_file <- file.path(root, "tables/BMC_Figure1_Main_MR_source.csv")
output_dir <- file.path(root, "figures")

dat <- read_csv(source_file, show_col_types = FALSE) |>
  mutate(
    Evidence_Class = case_when(
      Evidence_Tier == "Primary supported finding" ~ "Bonferroni-supported",
      Evidence_Tier == "Nominal finding" ~ "Nominal",
      grepl("Borderline", Evidence_Tier) ~ "Borderline / sensitivity-limited",
      TRUE ~ "Not statistically significant"
    ),
    Evidence_Class = factor(
      Evidence_Class,
      levels = c("Bonferroni-supported", "Nominal", "Borderline / sensitivity-limited", "Not statistically significant")
    )
  )

cols <- c(
  "Bonferroni-supported" = "#16679A",
  "Nominal" = "#2A9D8F",
  "Borderline / sensitivity-limited" = "#D08A00",
  "Not statistically significant" = "#7A7A7A"
)

fmt_p <- function(x) {
  sprintf("%.4f", x)
}

fmt_effect <- function(x, lo, hi, measure) {
  if (measure == "OR") {
    sprintf("%.3f (%.3f, %.3f)", x, lo, hi)
  } else if (measure == "Risk difference") {
    sprintf("%.4f (%.4f, %.5f)", x, lo, hi)
  } else {
    sprintf("%.3f (%.3f, %.3f)", x, lo, hi)
  }
}

text_theme <- theme_void(base_family = "Arial") +
  theme(
    plot.margin = margin(1, 2, 1, 2),
    plot.title = element_text(face = "bold", size = 9.3, hjust = 0, margin = margin(b = 2))
  )

panel_row <- function(d, log_scale = FALSE) {
  d <- d |>
    arrange(Outcome_Order) |>
    mutate(
      Outcome = factor(Outcome, levels = rev(Outcome)),
      effect_text = mapply(fmt_effect, Estimate, CI_Lower, CI_Upper, Effect_Measure),
      p_text = fmt_p(P)
    )

  label_plot <- ggplot(d, aes(y = Outcome, x = 0)) +
    geom_text(aes(label = as.character(Outcome)), hjust = 0, size = 3.15, family = "Arial") +
    xlim(0, 1) +
    text_theme

  effect_plot <- ggplot(d, aes(y = Outcome, x = 0)) +
    geom_text(aes(label = effect_text), hjust = 0, size = 3.0, family = "Arial") +
    xlim(0, 1) +
    text_theme

  forest <- ggplot(d, aes(y = Outcome, x = Estimate, xmin = CI_Lower, xmax = CI_Upper, colour = Evidence_Class)) +
    geom_vline(xintercept = ifelse(log_scale, 1, 0), linetype = "22", linewidth = 0.45, colour = "#555555") +
    geom_errorbarh(height = 0, linewidth = 0.72) +
    geom_point(size = 2.35) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = NULL, colour = NULL) +
    theme_classic(base_family = "Arial", base_size = 8.5) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(1, 3, 1, 3),
      legend.position = "none",
      legend.text = element_text(size = 7.8),
      legend.key.width = unit(9, "pt"),
      legend.spacing.x = unit(3, "pt")
    )

  if (log_scale) {
    forest <- forest +
      scale_x_log10(
        limits = c(0.12, 13),
        breaks = c(0.2, 0.5, 1, 2, 5, 10),
        labels = label_number(accuracy = 0.1)
      )
  } else if (all(d$Effect_Measure == "Risk difference")) {
    forest <- forest + scale_x_continuous(limits = c(-0.065, 0.015), breaks = c(-0.06, -0.03, 0))
  } else {
    forest <- forest + scale_x_continuous(limits = c(-0.20, 0.31), breaks = c(-0.2, 0, 0.2))
  }

  p_plot <- ggplot(d, aes(y = Outcome, x = 0)) +
    geom_text(aes(label = p_text), hjust = 0, size = 3.0, family = "Arial") +
    xlim(0, 1) +
    text_theme

  (label_plot | effect_plot | forest | p_plot) +
    plot_layout(widths = c(2.45, 2.25, 2.9, 0.9), guides = "collect")
}

column_header <- function(label) {
  ggplot() +
    annotate("text", x = 0, y = 0.5, label = label, hjust = 0, fontface = "bold", size = 3.35, family = "Arial") +
    xlim(0, 1) + ylim(0, 1) + theme_void()
}

section_header <- function(label) {
  ggplot() +
    annotate("text", x = 0, y = 0.5, label = label, hjust = 0, fontface = "bold", size = 3.55, family = "Arial") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(1, 2, 0, 2))
}

header <- (column_header("Outcome") | column_header("Estimate (95% CI)") |
  column_header("Forest plot") | column_header("P value")) +
  plot_layout(widths = c(2.45, 2.25, 2.9, 0.9))

p_a <- panel_row(filter(dat, Panel == "A_UKB_risk_difference"))
p_b <- panel_row(filter(dat, Panel == "B_FinnGen_OR"), log_scale = TRUE)
p_c <- panel_row(filter(dat, Panel == "C_lung_function_beta"))

legend_df <- tibble(
  x = c(0.8, 2.2, 3.25, 6.05),
  label = levels(dat$Evidence_Class),
  colour = unname(cols[levels(dat$Evidence_Class)])
)
legend_plot <- ggplot(legend_df, aes(x = x, y = 0.5)) +
  geom_point(aes(colour = label), size = 2.25) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.08, size = 2.75, family = "Arial") +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0.6, 8.3)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 5, 0, 5))

fig <- header /
  section_header("A  UK Biobank COPD differential diagnosis — risk difference") /
  p_a /
  section_header("B  FinnGen COPD outcomes — odds ratio") /
  p_b /
  section_header("C  Lung-function outcomes — beta coefficient") /
  p_c /
  legend_plot +
  plot_layout(heights = c(0.28, 0.30, 0.82, 0.30, 3.45, 0.30, 1.85, 0.42))

base <- file.path(output_dir, "Figure1_Main_Univariable_MR_FINAL")
width_in <- 183 / 25.4
height_in <- 154 / 25.4

svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(fig)
dev.off()

cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(fig)
dev.off()

agg_tiff(paste0(base, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
print(fig)
dev.off()

agg_png(paste0(base, ".png"), width = width_in, height = height_in, units = "in", res = 220)
print(fig)
dev.off()

