#!/usr/bin/env Rscript

# Every tree-building tool-pair correlation value, for the paper's concordance
# claim.
#
# The claim is about how much the choice of reconstruction tool changes the
# result, so the figure plots the individual pair values rather than a summary
# statistic. Two groupings are shown because they answer different questions and
# give different numbers:
#
#   Pooled       one correlation per tool pair over all 125 conditions, with the
#                five recorder technologies mixed together. Technologies differ
#                enormously in difficulty, so this measures largely whether both
#                tools agree that mitochondrial conditions are hard and BASELINE
#                conditions are easy. That shared signal inflates agreement.
#
#   Within       the same correlation computed separately inside each technology
#                (25 conditions per pair). Removing the between-technology
#                variance leaves only whether the tools rank the same recordings
#                the same way at a fixed recorder, which is the stricter and more
#                relevant question.
#
# Reported alongside is the same summary each grouping supports, so a caption can
# quote the number that matches its own qualifier.
#
# Usage:
#   Rscript plot-tool-concordance-values.R [results-dir] [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else {
  "test_harness/tree_building/results/recorder_parameter_grid_n250_level05"
}
root <- normalizePath(root, mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(root, "tool_concordance_values")
}

if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
  cache_dir <- file.path(tempdir(), "xdg-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = cache_dir)
}
for (package in c("ggplot2", "ragg")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ggplot2 and ragg packages are required to render the figure.")
  }
}

pooled_path <- file.path(
  root, "neutral_recorder_method_correlation_tool_correlation.tsv"
)
within_path <- file.path(
  root, "neutral_recorder_method_correlation_tool_correlation_by_system.tsv"
)
if (any(!file.exists(c(pooled_path, within_path)))) {
  stop(
    "Run correlate-neutral-recorder-methods.R first; missing:\n  ",
    paste(c(pooled_path, within_path)[!file.exists(c(pooled_path, within_path))],
          collapse = "\n  ")
  )
}
pooled <- read.delim(pooled_path, check.names = FALSE)
within <- read.delim(within_path, check.names = FALSE)

#' Keep one row per unordered tool pair
#'
#' The correlation tables tile a full square matrix, so each pair appears twice
#' and the diagonal is a self-correlation of 1. Both would distort a median.
#'
#' @param value Correlation rows carrying `method_1` and `method_2`.
#' @param extra Additional columns that form part of the grouping key.
#' @return The input restricted to unique unordered off-diagonal pairs, with a
#'   `pair_label` naming each one.
unique_pairs <- function(value, extra = character()) {
  value <- value[value$method_1 != value$method_2, , drop = FALSE]
  key <- paste(
    pmin(value$method_1, value$method_2),
    pmax(value$method_1, value$method_2),
    sep = "\r"
  )
  if (length(extra)) {
    key <- paste(key, do.call(paste, c(value[extra], sep = "\r")), sep = "\r")
  }
  value <- value[!duplicated(key), , drop = FALSE]
  value$pair_label <- paste(
    method_labels[pmin(value$method_1, value$method_2)],
    method_labels[pmax(value$method_1, value$method_2)],
    sep = " vs "
  )
  value
}

method_labels <- c(
  nj = "NJ", parsimony = "Parsimony", iqtree = "IQ-TREE 2",
  mix = "MIX", vine = "VINE", cassiopeia = "Cassiopeia"
)
recorder_labels <- c(
  baseline = "BASELINE", wt_crispr = "WT CRISPR (FLARE-like)",
  prime = "Prime editing", palincode = "PALINCODE",
  mitochondrial = "Mitochondrial"
)
recorder_levels <- names(recorder_labels)

pooled_pairs <- unique_pairs(pooled)
within_pairs <- unique_pairs(within, extra = "system")
within_pairs$system <- factor(
  within_pairs$system, levels = recorder_levels, labels = recorder_labels
)

pooled_median <- median(pooled_pairs$spearman_rho)
within_median <- median(within_pairs$spearman_rho)
# Order tool pairs by their within-technology median so the weakest pairs are
# immediately visible rather than buried in alphabetical order.
pair_order <- names(sort(tapply(
  within_pairs$spearman_rho, within_pairs$pair_label, median
)))
pooled_pairs$pair_label <- factor(pooled_pairs$pair_label, levels = pair_order)
within_pairs$pair_label <- factor(within_pairs$pair_label, levels = pair_order)

summary_table <- rbind(
  data.frame(
    grouping = "pooled_across_technologies",
    conditions_per_pair = max(pooled_pairs$n),
    pairs = nrow(pooled_pairs),
    median_rho = pooled_median,
    min_rho = min(pooled_pairs$spearman_rho),
    max_rho = max(pooled_pairs$spearman_rho),
    negative_pairs = sum(pooled_pairs$spearman_rho < 0),
    stringsAsFactors = FALSE
  ),
  data.frame(
    grouping = "within_each_technology",
    conditions_per_pair = max(within_pairs$n),
    pairs = nrow(within_pairs),
    median_rho = within_median,
    min_rho = min(within_pairs$spearman_rho),
    max_rho = max(within_pairs$spearman_rho),
    negative_pairs = sum(within_pairs$spearman_rho < 0),
    stringsAsFactors = FALSE
  )
)
summary_path <- paste0(output_prefix, "_summary.tsv")
write.table(
  summary_table, summary_path, sep = "\t", row.names = FALSE, quote = FALSE
)
values_path <- paste0(output_prefix, "_values.tsv")
write.table(
  within_pairs[c(
    "system", "method_1", "method_2", "n", "spearman_rho", "pearson_r"
  )],
  values_path, sep = "\t", row.names = FALSE, quote = FALSE
)

ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 12, base_family = "sans")
)
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 13),
  plot.subtitle = ggplot2::element_text(size = 9.5, color = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank()
)

# Panel 1. Every value, one row per tool pair, coloured by technology, with the
# pooled value marked separately so the two groupings can be read against each
# other on one axis.
combined_points <- rbind(
  data.frame(
    pair_label = within_pairs$pair_label, grouping = as.character(within_pairs$system),
    spearman_rho = within_pairs$spearman_rho, stringsAsFactors = FALSE
  ),
  data.frame(
    pair_label = pooled_pairs$pair_label, grouping = "Pooled over all technologies",
    spearman_rho = pooled_pairs$spearman_rho, stringsAsFactors = FALSE
  )
)
grouping_levels <- c(unname(recorder_labels), "Pooled over all technologies")
combined_points$grouping <- factor(combined_points$grouping, levels = grouping_levels)
grouping_colors <- c(
  BASELINE = "#B4436C", `WT CRISPR (FLARE-like)` = "#DC6B2F",
  `Prime editing` = "#087E8B", PALINCODE = "#708B36",
  Mitochondrial = "#3D67A8", `Pooled over all technologies` = "#1F3B4D"
)
grouping_shapes <- c(
  BASELINE = 16, `WT CRISPR (FLARE-like)` = 16, `Prime editing` = 16,
  PALINCODE = 16, Mitochondrial = 16, `Pooled over all technologies` = 18
)

value_panel <- ggplot2::ggplot(
  combined_points, ggplot2::aes(x = spearman_rho, y = pair_label)
) +
  ggplot2::geom_vline(xintercept = 0, color = "#8A8A8A", linewidth = 0.4) +
  ggplot2::geom_vline(
    xintercept = pooled_median, color = "#1F3B4D",
    linewidth = 0.5, linetype = "dashed"
  ) +
  ggplot2::geom_vline(
    xintercept = within_median, color = "#B4436C",
    linewidth = 0.5, linetype = "dashed"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = grouping, shape = grouping, size = grouping),
    alpha = 0.9
  ) +
  ggplot2::scale_color_manual(name = NULL, values = grouping_colors) +
  ggplot2::scale_shape_manual(name = NULL, values = grouping_shapes) +
  ggplot2::scale_size_manual(
    name = NULL,
    values = c(rep(2.3, length(recorder_labels)), 3.4)
  ) +
  ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  ggplot2::labs(
    title = "Every tool-pair correlation, by grouping",
    subtitle = sprintf(
      paste(
        "Diamonds pool all technologies (dashed navy, median %.2f);",
        "circles are within one technology (dashed pink, median %.2f)."
      ),
      pooled_median, within_median
    ),
    x = "Spearman rho of normalized RF", y = NULL
  ) +
  panel_theme +
  ggplot2::theme(legend.position = "right")

# Panel 2. The two groupings as distributions, which is the comparison the
# caption's number depends on.
distribution_points <- rbind(
  data.frame(grouping = "Within one\ntechnology",
             spearman_rho = within_pairs$spearman_rho, stringsAsFactors = FALSE),
  data.frame(grouping = "Pooled over all\ntechnologies",
             spearman_rho = pooled_pairs$spearman_rho, stringsAsFactors = FALSE)
)
distribution_points$grouping <- factor(
  distribution_points$grouping,
  levels = c("Pooled over all\ntechnologies", "Within one\ntechnology")
)
distribution_panel <- ggplot2::ggplot(
  distribution_points, ggplot2::aes(x = grouping, y = spearman_rho)
) +
  ggplot2::geom_hline(yintercept = 0, color = "#8A8A8A", linewidth = 0.4) +
  ggplot2::geom_violin(
    ggplot2::aes(fill = grouping), color = NA, alpha = 0.35, scale = "width"
  ) +
  ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "#FFFFFF") +
  ggplot2::geom_jitter(
    width = 0.07, height = 0, size = 1.5, alpha = 0.65, color = "#333333"
  ) +
  ggplot2::scale_fill_manual(
    values = c(`Pooled over all\ntechnologies` = "#1F3B4D",
               `Within one\ntechnology` = "#B4436C"),
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  ggplot2::labs(
    title = "The grouping changes the number",
    subtitle = "Each point is one tool pair",
    x = NULL, y = "Spearman rho"
  ) +
  panel_theme

# Panel 3. Within-technology values split by recorder, which locates the
# disagreement rather than averaging over it.
system_panel <- ggplot2::ggplot(
  within_pairs, ggplot2::aes(x = system, y = spearman_rho)
) +
  ggplot2::geom_hline(yintercept = 0, color = "#8A8A8A", linewidth = 0.4) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = system), width = 0.55, alpha = 0.35,
    color = "#333333", outlier.shape = NA
  ) +
  ggplot2::geom_jitter(width = 0.09, height = 0, size = 1.6, alpha = 0.75) +
  ggplot2::scale_fill_manual(values = unname(grouping_colors[1:5]), guide = "none") +
  ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  ggplot2::labs(
    title = "Within-technology agreement, by recorder",
    subtitle = "The negative values are Cassiopeia against every other tool on BASELINE",
    x = NULL, y = "Spearman rho"
  ) +
  panel_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(
    value_panel,
    patchwork::wrap_plots(distribution_panel, system_panel, ncol = 2,
                          widths = c(0.75, 1)),
    ncol = 1, heights = c(1, 0.85)
  ) +
    patchwork::plot_annotation(
      title = "Concordance between tree-building tools on neutral trees",
      subtitle = sprintf(
        "%d tool pairs; %d conditions per pooled pair, %d per within-technology pair",
        nrow(pooled_pairs), max(pooled_pairs$n), max(within_pairs$n)
      ),
      # Panel tags so a figure legend can refer to (A), (B), (C).
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#4A4A4A"),
        plot.tag = ggplot2::element_text(face = "bold", size = 15)
      )
    )
} else {
  value_panel
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(
  png_path, plot = combined, device = ragg::agg_png,
  width = 12.5, height = 10, units = "in", dpi = 220, background = "#FFFFFF"
)
ggplot2::ggsave(
  pdf_path, plot = combined, device = grDevices::cairo_pdf,
  width = 12.5, height = 10, units = "in", bg = "#FFFFFF"
)

print(summary_table, row.names = FALSE)
cat("\nWithin-technology median by recorder:\n")
by_system <- tapply(within_pairs$spearman_rho, within_pairs$system, median)
for (name in names(by_system)) {
  cat(sprintf("  %-24s %.3f\n", name, by_system[[name]]))
}
negative <- within_pairs[within_pairs$spearman_rho < 0, , drop = FALSE]
if (nrow(negative)) {
  cat("\nNegative within-technology pairs:\n")
  print(
    negative[order(negative$spearman_rho),
             c("system", "method_1", "method_2", "n", "spearman_rho")],
    row.names = FALSE
  )
}
cat("\nWrote:\n", paste0("  ", c(png_path, pdf_path, summary_path, values_path),
                         collapse = "\n"), "\n")
