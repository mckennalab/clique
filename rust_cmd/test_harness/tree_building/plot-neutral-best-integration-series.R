#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else {
  "test_harness/tree_building/results/neutral_best_integration_series_n250"
}
root <- normalizePath(root, mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(root, "neutral_best_integration_series")
}

if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
  cache_dir <- file.path(tempdir(), "xdg-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = cache_dir)
}
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("ragg", quietly = TRUE)) {
  stop("The ggplot2 and ragg packages are required to render the figure.")
}

overall_path <- file.path(root, "accuracy_by_recorder_integration.tsv")
method_path <- file.path(root, "accuracy_by_recorder_integration_method.tsv")
if (!file.exists(overall_path) || !file.exists(method_path)) {
  stop("Run the neutral integration series before plotting.")
}
overall <- read.delim(overall_path, check.names = FALSE)
methods <- read.delim(method_path, check.names = FALSE)
if (!nrow(overall) || !nrow(methods)) {
  stop("The neutral integration summaries are empty.")
}

recorder_levels <- c(
  "baseline", "wt_crispr", "prime", "palincode", "mitochondrial"
)
recorder_labels <- c(
  baseline = "BASELINE\nbest rate: 0.00952",
  wt_crispr = "WT CRISPR (FLARE-like)\nbest rate: 0.005 / 0.05 ins/del",
  prime = "Prime editing\nbest rate: 0.04",
  palincode = "PALINCODE\nbest rate: 0.08 / 0.05",
  mitochondrial = "Mitochondrial\nbest rate: 8e-6"
)
method_levels <- c("nj", "parsimony", "iqtree", "mix", "vine", "cassiopeia")
method_labels <- c(
  nj = "NJ",
  parsimony = "Basic parsimony",
  iqtree = "IQ-TREE 2",
  mix = "MIX",
  vine = "VINE",
  cassiopeia = "Cassiopeia greedy"
)
method_colors <- c(
  nj = "#087E8B",
  parsimony = "#DC6B2F",
  iqtree = "#3D67A8",
  mix = "#D29A18",
  vine = "#708B36",
  cassiopeia = "#9B5575"
)

overall$system <- factor(overall$system, levels = recorder_levels)
methods$system <- factor(methods$system, levels = recorder_levels)
methods$method <- factor(methods$method, levels = method_levels)
if (anyNA(overall$system) || anyNA(methods$system) || anyNA(methods$method)) {
  stop("The summaries contain an unexpected recorder or reconstruction method.")
}
thresholds <- unique(overall$threshold)
if (length(thresholds) != 1L || !is.finite(thresholds)) {
  stop("The summaries must contain one finite stopping threshold.")
}
threshold <- thresholds[[1L]]
level_breaks <- sort(unique(c(overall$level, methods$level)))

plot_value <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = threshold,
    color = "#A54232",
    linetype = "22",
    linewidth = 0.8
  ) +
  ggplot2::geom_errorbar(
    data = methods,
    ggplot2::aes(
      x = level,
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf,
      color = method,
      group = method
    ),
    width = 0.55,
    linewidth = 0.48,
    alpha = 0.62
  ) +
  ggplot2::geom_line(
    data = methods,
    ggplot2::aes(
      x = level,
      y = mean_normalized_rf,
      color = method,
      group = method
    ),
    linewidth = 0.72,
    alpha = 0.80
  ) +
  ggplot2::geom_point(
    data = methods,
    ggplot2::aes(x = level, y = mean_normalized_rf, color = method),
    size = 2.0,
    alpha = 0.92
  ) +
  ggplot2::geom_errorbar(
    data = overall,
    ggplot2::aes(
      x = level,
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf
    ),
    width = 0.9,
    linewidth = 0.82,
    color = "#171717"
  ) +
  ggplot2::geom_line(
    data = overall,
    ggplot2::aes(x = level, y = mean_normalized_rf, group = system),
    color = "#171717",
    linewidth = 1.25
  ) +
  ggplot2::geom_point(
    data = overall,
    ggplot2::aes(x = level, y = mean_normalized_rf),
    shape = 21,
    size = 3.5,
    stroke = 1.0,
    color = "#171717",
    fill = "#FFFDF7"
  ) +
  ggplot2::geom_point(
    data = overall[overall$threshold_reached, , drop = FALSE],
    ggplot2::aes(x = level, y = mean_normalized_rf),
    shape = 23,
    size = 4.8,
    stroke = 1.0,
    color = "#171717",
    fill = "#E9C75C"
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(system),
    ncol = 2,
    labeller = ggplot2::as_labeller(recorder_labels)
  ) +
  ggplot2::scale_x_continuous(
    breaks = level_breaks,
    minor_breaks = NULL,
    expand = ggplot2::expansion(mult = c(0.04, 0.06))
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    limits = c(0, 1.02),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::scale_color_manual(
    values = method_colors,
    breaks = method_levels,
    labels = method_labels
  ) +
  ggplot2::labs(
    title = "Integration sensitivity at recorder-specific best editing rates",
    subtitle = paste(
      "Neutral lineage trees with 250 sampled cells;",
      sprintf(
        "each recorder stops after mean normalized RF falls below %.0f%%",
        100 * threshold
      )
    ),
    x = "Integrations (recovered genomes for mitochondrial recording)",
    y = "Mean normalized RF distance",
    color = "Tree-building method",
    caption = paste(
      "Colored lines: method means. Black line: mean across all six methods.",
      sprintf(
        "Error bars: 95%% CIs across population seeds. Red dashed line: %.0f%% threshold.",
        100 * threshold
      ),
      "Gold diamond: stopping point."
    )
  ) +
  ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "#FFFFFF", color = NA),
    panel.background = ggplot2::element_rect(fill = "#FFFFFF", color = NA),
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      fill = "#FFFFFF", color = "#4F5D5F", linewidth = 0.7
    ),
    strip.text = ggplot2::element_text(
      color = "#384849", face = "bold", size = 10.5, lineheight = 1.1,
      margin = ggplot2::margin(8, 6, 8, 6)
    ),
    plot.title = ggplot2::element_text(
      color = "#173A40", face = "bold", size = 20
    ),
    plot.subtitle = ggplot2::element_text(color = "#4F5D5F", size = 11.5),
    plot.caption = ggplot2::element_text(
      color = "#5D6260", size = 9, hjust = 0, lineheight = 1.1
    ),
    axis.title = ggplot2::element_text(color = "#263536", face = "bold"),
    axis.text = ggplot2::element_text(color = "#384849"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    legend.key.width = grid::unit(1.3, "lines"),
    panel.spacing = grid::unit(0.9, "lines"),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  )

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  png_path,
  plot = plot_value,
  device = ragg::agg_png,
  width = 12,
  height = 10.5,
  units = "in",
  dpi = 300,
  background = "#FFFFFF"
)
ggplot2::ggsave(
  pdf_path,
  plot = plot_value,
  device = grDevices::cairo_pdf,
  width = 12,
  height = 10.5,
  units = "in",
  bg = "#FFFFFF"
)
cat("PNG: ", normalizePath(png_path), "\n", sep = "")
cat("PDF: ", normalizePath(pdf_path), "\n", sep = "")
