#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else {
  "test_harness/tree_building/results/neutral_best_nj_integration_series_nested_n250"
}
root <- normalizePath(root, mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(root, "neutral_best_nj_integration_series")
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

accuracy_path <- file.path(root, "accuracy_by_recorder_integration.tsv")
if (!file.exists(accuracy_path)) {
  stop("Run the NJ-only neutral integration series before plotting.")
}
accuracy <- read.delim(accuracy_path, check.names = FALSE)
if (!nrow(accuracy) || any(accuracy$best_method != "nj")) {
  stop("The accuracy table does not contain a complete NJ-only series.")
}

recorder_levels <- c(
  "baseline", "wt_crispr", "prime", "palincode", "mitochondrial"
)
recorder_labels <- c(
  baseline = "BASELINE\nbest rate: 0.00952",
  wt_crispr = "WT CRISPR (FLARE-like)\nbest rate: 0.005 / 0.05 ins/del",
  # PEtracer: 3 edit sites per integration, each installing one of 8 marks.
  # Kept in step with best_settings in generate-recorder-parameter-grid.R and
  # parameter_labels in finalize-neutral-best-integration-series.R.
  prime = "Prime editing (PEtracer)\nbest rate: 0.18, 3 sites x 8 marks",
  palincode = "PALINCODE\nbest rate: 0.08 / 0.05",
  mitochondrial = "Mitochondrial\nbest rate: 8e-6"
)
accuracy$system <- factor(accuracy$system, levels = recorder_levels)
if (anyNA(accuracy$system)) {
  stop("The accuracy table contains an unexpected recorder.")
}
thresholds <- unique(accuracy$threshold)
if (length(thresholds) != 1L || !is.finite(thresholds)) {
  stop("The accuracy table must contain one finite threshold.")
}
threshold <- thresholds[[1L]]
level_breaks <- sort(unique(accuracy$level))

reached_rows <- accuracy[accuracy$threshold_reached, , drop = FALSE]
stopping_rows <- if (nrow(reached_rows)) {
  do.call(rbind, lapply(
    split(reached_rows, reached_rows$system, drop = TRUE),
    function(part) part[which.min(part$level), , drop = FALSE]
  ))
} else {
  accuracy[FALSE, , drop = FALSE]
}

plot_value <- ggplot2::ggplot(
  accuracy,
  ggplot2::aes(x = level, y = mean_normalized_rf, group = system)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf
    ),
    width = 2.0,
    linewidth = 0.75,
    color = "#087E8B"
  ) +
  ggplot2::geom_line(color = "#087E8B", linewidth = 1.25) +
  ggplot2::geom_point(
    shape = 21,
    size = 4.0,
    stroke = 1.0,
    color = "#075965",
    fill = "#FFFDF7"
  ) +
  ggplot2::geom_point(
    data = stopping_rows,
    shape = 23,
    size = 5.2,
    stroke = 1.0,
    color = "#171717",
    fill = "#E9C75C"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.3f", mean_normalized_rf)),
    color = "#263536",
    fontface = "bold",
    size = 3.2,
    vjust = -1.25
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(system),
    ncol = 2,
    labeller = ggplot2::as_labeller(recorder_labels)
  ) +
  ggplot2::scale_x_continuous(
    breaks = level_breaks,
    minor_breaks = NULL,
    expand = ggplot2::expansion(mult = c(0.08, 0.10))
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    limits = c(0, 1.08),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "NJ accuracy with increasing recorder integrations",
    subtitle = paste(
      "Nested subsets of one maximum-level recording per replicate;",
      "neutral trees with 250 sampled cells;",
      "100 mitochondrial genomes simulated per cell"
    ),
    x = "Integrations (recovered genomes for mitochondrial recording)",
    y = "Mean normalized RF distance",
    caption = sprintf(
      paste(
        "Points: NJ means across five population seeds. Error bars: 95%% CIs.",
        "Gold diamond: first requested level below the %.0f%% target."
      ),
      100 * threshold
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
    panel.spacing = grid::unit(0.9, "lines"),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
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
