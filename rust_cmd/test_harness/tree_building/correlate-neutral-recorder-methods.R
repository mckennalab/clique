#!/usr/bin/env Rscript

# Concordance between tree-building tools and between recorder technologies on
# one tree shape (neutral by default).
#
# The recorder-parameter grid reconstructs every condition with every tool, so
# the same condition appears once per tool and the same (tool, tier, seed)
# appears once per technology. That pairing supports two questions:
#
#   Tool concordance      Do two tools rank the same conditions as easy and
#                         hard? Pairs are matched on condition_id, so each
#                         correlation is over the shared set of conditions and
#                         asks whether tool choice reorders the results. A high
#                         value means the tools agree about which recordings
#                         carry lineage signal, even if their absolute accuracy
#                         differs; a low value means a conclusion drawn with one
#                         tool would not survive swapping the tool.
#
#   Technology concordance
#                         Do two recorders succeed and fail on the same
#                         (tool, editing tier, population seed)? Pairs are
#                         matched on that triple. Technologies do not share an
#                         editing rate, so they are aligned by parameter_rank
#                         (tier position), not by rate value.
#
# Correlations are computed on normalized RF. Spearman is the headline statistic
# because the question is about ordering and the accuracy distributions are
# skewed and bounded; Pearson is reported alongside it. Both are identical
# whether computed on normalized RF or on accuracy (1 - normalized RF), since
# that transform flips both series.
#
# Correlation here measures agreement in *ordering*, not agreement in value. Two
# tools can correlate at 1.0 while one is uniformly far less accurate, so read
# these next to the accuracy tables rather than instead of them.
#
# Usage:
#   Rscript correlate-neutral-recorder-methods.R [results-root] [shape] [output-prefix]
#
# Writes <output-prefix>.{png,pdf} and, next to them, three TSVs:
#   *_tool_correlation.tsv             tool x tool, pooled over technologies
#   *_tool_correlation_by_system.tsv   tool x tool, within each technology
#   *_technology_correlation.tsv       technology x technology

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else {
  "test_harness/tree_building/results/recorder_parameter_grid_n250_level05"
}
root <- normalizePath(root, mustWork = TRUE)
shape_filter <- if (length(arguments) >= 2L) tolower(arguments[[2L]]) else "neutral"
output_prefix <- if (length(arguments) >= 3L) arguments[[3L]] else {
  file.path(root, sprintf("%s_recorder_method_correlation", shape_filter))
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

results_path <- file.path(root, "reconstruction_with_parameters.tsv.gz")
if (!file.exists(results_path)) {
  stop("Run summarize-recorder-parameter-grid.R before correlating results.")
}
results <- read.delim(gzfile(results_path), check.names = FALSE)
required_columns <- c(
  "condition_id", "shape", "population_seed", "system", "method",
  "parameter_rank", "normalized_rf", "status"
)
if (any(!required_columns %in% names(results))) {
  stop("reconstruction_with_parameters.tsv.gz is missing required columns.")
}

results <- results[tolower(results$shape) == shape_filter, , drop = FALSE]
if (!nrow(results)) {
  stop(sprintf("No reconstructions are present for shape %s.", shape_filter))
}
if (any(results$status != "success")) {
  stop("Refusing to correlate a grid that contains failed reconstructions.")
}
if (anyNA(results$normalized_rf)) {
  stop("Refusing to correlate a grid that contains missing normalized RF values.")
}

method_levels <- c("nj", "parsimony", "iqtree", "mix", "vine", "cassiopeia")
method_labels <- c(
  nj = "NJ", parsimony = "Basic\nparsimony", iqtree = "IQ-TREE 2",
  mix = "MIX", vine = "VINE", cassiopeia = "Cassiopeia\ngreedy"
)
recorder_levels <- c(
  "baseline", "wt_crispr", "prime", "palincode", "mitochondrial"
)
recorder_labels <- c(
  baseline = "BASELINE", wt_crispr = "WT CRISPR\n(FLARE-like)",
  prime = "Prime\nediting", palincode = "PALINCODE",
  mitochondrial = "Mito-\nchondrial"
)
# Single-line forms for facet strips, where a wrapped label reads as two words.
recorder_flat_labels <- c(
  baseline = "BASELINE", wt_crispr = "WT CRISPR (FLARE-like)",
  prime = "Prime editing", palincode = "PALINCODE",
  mitochondrial = "Mitochondrial"
)
method_levels <- method_levels[method_levels %in% results$method]
recorder_levels <- recorder_levels[recorder_levels %in% results$system]
if (length(method_levels) < 2L || length(recorder_levels) < 2L) {
  stop("At least two tools and two technologies are required to correlate.")
}

#' Correlate one measured column between every pair of groups
#'
#' Reshapes long results into a group-by-key matrix and correlates each pair of
#' groups over the keys they share. Pairs are matched on `key_columns`, so a
#' group missing a key contributes nothing for that key rather than shifting the
#' alignment.
#'
#' @param value Long data frame of reconstruction results.
#' @param group_column Column naming the series to correlate (a tool or a
#'   technology).
#' @param key_columns Columns whose combination identifies one paired
#'   observation shared across groups.
#' @param group_levels Groups to include, in output order.
#' @return Data frame with one row per ordered group pair: both group names,
#'   the number of shared keys `n`, `spearman_rho`, `pearson_r`, and the
#'   two-sided Spearman p-value `spearman_p`. Self-pairs are included with a
#'   correlation of 1 so the result tiles a complete square matrix.
correlate_groups <- function(value, group_column, key_columns, group_levels) {
  value$pair_key <- do.call(paste, c(value[key_columns], sep = "\r"))
  series <- lapply(group_levels, function(group) {
    part <- value[value[[group_column]] == group, , drop = FALSE]
    if (anyDuplicated(part$pair_key)) {
      stop(sprintf(
        "Group %s has more than one result per pairing key; the key columns do not identify a unique observation.",
        group
      ))
    }
    stats::setNames(part$normalized_rf, part$pair_key)
  })
  names(series) <- group_levels

  rows <- list()
  for (first in group_levels) {
    for (second in group_levels) {
      left <- series[[first]]
      right <- series[[second]]
      shared <- intersect(names(left), names(right))
      left_values <- unname(left[shared])
      right_values <- unname(right[shared])
      # A constant series has no ranking to agree with, so correlation is
      # undefined rather than zero; leave it NA and let the plot show a gap.
      degenerate <- length(shared) < 3L ||
        stats::sd(left_values) == 0 || stats::sd(right_values) == 0
      spearman <- if (identical(first, second)) {
        1
      } else if (degenerate) {
        NA_real_
      } else {
        suppressWarnings(stats::cor(
          left_values, right_values, method = "spearman"
        ))
      }
      pearson <- if (identical(first, second)) {
        1
      } else if (degenerate) {
        NA_real_
      } else {
        stats::cor(left_values, right_values, method = "pearson")
      }
      spearman_p <- if (identical(first, second) || degenerate) {
        NA_real_
      } else {
        suppressWarnings(stats::cor.test(
          left_values, right_values, method = "spearman"
        )$p.value)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        first = first, second = second, n = length(shared),
        spearman_rho = spearman, pearson_r = pearson,
        spearman_p = spearman_p, stringsAsFactors = FALSE
      )
    }
  }
  output <- do.call(rbind, rows)
  names(output)[1:2] <- paste0(group_column, c("_1", "_2"))
  output
}

# Tools are paired on the condition they reconstructed: same recording matrix,
# different algorithm.
tool_correlation <- correlate_groups(
  results, "method", "condition_id", method_levels
)
tool_correlation$shape <- shape_filter

# The same comparison inside each technology, which separates "these tools
# agree" from "these tools agree only because technologies differ so much".
tool_correlation_by_system <- do.call(rbind, lapply(
  recorder_levels,
  function(system_name) {
    part <- results[results$system == system_name, , drop = FALSE]
    rows <- correlate_groups(part, "method", "condition_id", method_levels)
    rows$system <- system_name
    rows$shape <- shape_filter
    rows
  }
))

# Technologies are paired on tool, editing tier, and population seed. Tier
# position aligns them because the technologies do not share a rate scale.
technology_correlation <- correlate_groups(
  results, "system",
  c("method", "parameter_rank", "population_seed"),
  recorder_levels
)
technology_correlation$shape <- shape_filter

tool_path <- paste0(output_prefix, "_tool_correlation.tsv")
tool_system_path <- paste0(output_prefix, "_tool_correlation_by_system.tsv")
technology_path <- paste0(output_prefix, "_technology_correlation.tsv")
write_tsv <- function(value, path) {
  write.table(value, path, sep = "\t", row.names = FALSE, quote = FALSE)
}
write_tsv(tool_correlation, tool_path)
write_tsv(tool_correlation_by_system, tool_system_path)
write_tsv(technology_correlation, technology_path)

#' Build one annotated correlation heatmap
#'
#' @param value Correlation rows from `correlate_groups()`.
#' @param first_column,second_column Names of the two group columns.
#' @param level_labels Named vector of display labels, ordered as the axes
#'   should read.
#' @param title,subtitle Panel text.
#' @return A ggplot object.
correlation_heatmap <- function(value, first_column, second_column,
                                level_labels, title, subtitle) {
  level_names <- names(level_labels)
  value$x <- factor(value[[first_column]], levels = level_names)
  value$y <- factor(value[[second_column]], levels = rev(level_names))
  value$label <- ifelse(
    is.na(value$spearman_rho), "—", sprintf("%.2f", value$spearman_rho)
  )
  # Label contrast is chosen from the fill, not from the value's sign, so it
  # stays readable at both ends of the scale.
  value$label_color <- ifelse(
    !is.na(value$spearman_rho) & abs(value$spearman_rho) > 0.6,
    "#FFFFFF", "#1A1A1A"
  )
  ggplot2::ggplot(value, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = spearman_rho), color = "#FFFFFF", linewidth = 0.9
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, color = label_color),
      size = 3.4, fontface = "bold"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(
      name = "Spearman rho",
      low = "#3D67A8", mid = "#F4F4F2", high = "#B4436C",
      midpoint = 0, limits = c(-1, 1), na.value = "#E4E4E0"
    ) +
    ggplot2::scale_x_discrete(labels = level_labels, position = "top") +
    ggplot2::scale_y_discrete(labels = level_labels) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9, lineheight = 0.95),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.95),
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9.5, color = "#4A4A4A")
    )
}

shape_title <- paste0(
  toupper(substring(shape_filter, 1, 1)), substring(shape_filter, 2)
)
tool_panel <- correlation_heatmap(
  tool_correlation, "method_1", "method_2", method_labels[method_levels],
  sprintf("Tree-building tools agree on which %s recordings are hard", shape_filter),
  sprintf(
    "Paired on the %d shared conditions\n(%d technologies x %d tiers x %d seeds)",
    max(tool_correlation$n),
    length(recorder_levels),
    length(unique(results$parameter_rank)),
    length(unique(results$population_seed))
  )
)
technology_panel <- correlation_heatmap(
  technology_correlation, "system_1", "system_2",
  recorder_labels[recorder_levels],
  "Recorder technologies agree far less with each other",
  sprintf(
    "Paired on %d shared tool-by-tier-by-seed\ncombinations",
    max(technology_correlation$n)
  )
)

# Pooling tools over all technologies inflates their agreement: technologies
# differ so much in difficulty that every tool agrees mitochondrial conditions
# are hard and BASELINE conditions are easy, before any tool-specific behavior
# enters. Repeating the comparison inside each technology removes that shared
# variance and is the honest measure of whether tool choice matters.
within_system <- tool_correlation_by_system
within_system$x <- factor(within_system$method_1, levels = method_levels)
within_system$y <- factor(within_system$method_2, levels = rev(method_levels))
within_system$system <- factor(within_system$system, levels = recorder_levels)
within_system$label <- ifelse(
  is.na(within_system$spearman_rho), "—",
  sprintf("%.2f", within_system$spearman_rho)
)
within_panel <- ggplot2::ggplot(within_system, ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_tile(
    ggplot2::aes(fill = spearman_rho), color = "#FFFFFF", linewidth = 0.6
  ) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5) +
  ggplot2::scale_fill_gradient2(
    name = "Spearman rho", low = "#3D67A8", mid = "#F4F4F2",
    high = "#B4436C", midpoint = 0, limits = c(-1, 1), na.value = "#E4E4E0"
  ) +
  ggplot2::scale_x_discrete(
    labels = gsub("\n", " ", method_labels[method_levels])
  ) +
  ggplot2::scale_y_discrete(labels = gsub("\n", " ", method_labels[method_levels])) +
  ggplot2::facet_wrap(
    ggplot2::vars(system), nrow = 1,
    labeller = ggplot2::as_labeller(recorder_flat_labels[recorder_levels])
  ) +
  ggplot2::labs(
    title = "Within a single technology, tools agree much less",
    subtitle = paste(
      "One technology at a time, 25 conditions each (5 rate tiers x 5 seeds).",
      "Cassiopeia is anti-correlated with every other tool on BASELINE."
    ),
    x = NULL, y = NULL
  ) +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 6.5, angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 6.5),
    strip.text = ggplot2::element_text(face = "bold", size = 9),
    plot.title = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle = ggplot2::element_text(size = 9.5, color = "#4A4A4A")
  )

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  top_row <- patchwork::wrap_plots(
    tool_panel, technology_panel, ncol = 2, widths = c(1, 1)
  )
  patchwork::wrap_plots(top_row, within_panel, nrow = 2, heights = c(1, 0.85)) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf(
        "%s trees: concordance across tree-building tools and recorder technologies",
        shape_title
      ),
      caption = paste(
        "Spearman rho on normalized Robinson-Foulds distance; identical to",
        "using accuracy = 1 - normalized RF.\nCorrelation measures agreement in",
        "ordering, not in absolute accuracy."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15),
        plot.caption = ggplot2::element_text(size = 9, color = "#4A4A4A", hjust = 0)
      )
    )
} else {
  # patchwork is optional; fall back to a single stacked figure built from the
  # long tables so the script still produces a plot without it.
  tool_rows <- tool_correlation
  names(tool_rows)[1:2] <- c("group_1", "group_2")
  tool_rows$panel <- "Tree-building tools"
  technology_rows <- technology_correlation
  names(technology_rows)[1:2] <- c("group_1", "group_2")
  technology_rows$panel <- "Recorder technologies"
  stacked <- rbind(
    tool_rows[c("group_1", "group_2", "spearman_rho", "panel")],
    technology_rows[c("group_1", "group_2", "spearman_rho", "panel")]
  )
  all_labels <- c(method_labels[method_levels], recorder_labels[recorder_levels])
  stacked$x <- factor(stacked$group_1, levels = names(all_labels))
  stacked$y <- factor(stacked$group_2, levels = rev(names(all_labels)))
  stacked$label <- ifelse(
    is.na(stacked$spearman_rho), "—", sprintf("%.2f", stacked$spearman_rho)
  )
  ggplot2::ggplot(stacked, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = spearman_rho), color = "#FFFFFF", linewidth = 0.9
    ) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2) +
    ggplot2::scale_fill_gradient2(
      name = "Spearman rho", low = "#3D67A8", mid = "#F4F4F2",
      high = "#B4436C", midpoint = 0, limits = c(-1, 1), na.value = "#E4E4E0"
    ) +
    ggplot2::scale_x_discrete(labels = all_labels, position = "top") +
    ggplot2::scale_y_discrete(labels = all_labels) +
    ggplot2::facet_wrap(ggplot2::vars(panel), scales = "free") +
    ggplot2::labs(
      title = sprintf("%s trees: tool and technology concordance", shape_title),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(
  png_path, plot = combined, device = ragg::agg_png,
  width = 13.5, height = 11.5, units = "in", dpi = 300, background = "#FFFFFF"
)
ggplot2::ggsave(
  pdf_path, plot = combined, device = grDevices::cairo_pdf,
  width = 13.5, height = 11.5, units = "in", bg = "#FFFFFF"
)

off_diagonal <- function(value, first_column, second_column) {
  value[value[[first_column]] != value[[second_column]], , drop = FALSE]
}
tool_off <- off_diagonal(tool_correlation, "method_1", "method_2")
technology_off <- off_diagonal(technology_correlation, "system_1", "system_2")

cat(sprintf("Shape: %s\n", shape_filter))
cat(sprintf(
  "Tool concordance (Spearman rho over %d shared conditions): median %.3f, range %.3f to %.3f\n",
  max(tool_correlation$n),
  median(tool_off$spearman_rho, na.rm = TRUE),
  min(tool_off$spearman_rho, na.rm = TRUE),
  max(tool_off$spearman_rho, na.rm = TRUE)
))
cat(sprintf(
  "Technology concordance (Spearman rho over %d shared combinations): median %.3f, range %.3f to %.3f\n",
  max(technology_correlation$n),
  median(technology_off$spearman_rho, na.rm = TRUE),
  min(technology_off$spearman_rho, na.rm = TRUE),
  max(technology_off$spearman_rho, na.rm = TRUE)
))
cat("\nTool pairs, least concordant first:\n")
tool_print <- tool_off[order(tool_off$spearman_rho), , drop = FALSE]
print(
  head(tool_print[c("method_1", "method_2", "n", "spearman_rho", "pearson_r")], 6),
  row.names = FALSE
)
cat("\nWrote:\n")
cat(paste0("  ", c(png_path, pdf_path, tool_path, tool_system_path, technology_path), collapse = "\n"), "\n")
