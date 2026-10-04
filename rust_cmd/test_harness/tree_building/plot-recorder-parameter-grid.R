#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else {
  "test_harness/tree_building/results/recorder_parameter_grid_n250_level05"
}
root <- normalizePath(root, mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(root, "recorder_parameter_grid_summary")
}
shape_output_prefix <- if (length(arguments) >= 3L) arguments[[3L]] else {
  file.path(root, "recorder_parameter_grid_by_tree_shape")
}
neutral_output_prefix <- if (length(arguments) >= 4L) arguments[[4L]] else {
  file.path(root, "recorder_parameter_grid_neutral")
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

parameter_path <- file.path(root, "accuracy_by_recorder_parameter.tsv")
method_path <- file.path(root, "accuracy_by_recorder_parameter_method.tsv")
shape_path <- file.path(root, "accuracy_by_recorder_parameter_shape.tsv")
shape_method_path <- file.path(
  root, "accuracy_by_recorder_parameter_shape_method.tsv"
)
required_paths <- c(parameter_path, method_path, shape_path, shape_method_path)
if (any(!file.exists(required_paths))) {
  stop("Run summarize-recorder-parameter-grid.R before plotting.")
}

overall <- read.delim(parameter_path, check.names = FALSE)
methods <- read.delim(method_path, check.names = FALSE)
by_shape <- read.delim(shape_path, check.names = FALSE)
methods_by_shape <- read.delim(shape_method_path, check.names = FALSE)
expected_method_rows <- nrow(overall) * 6L
if (!nrow(overall) || nrow(methods) != expected_method_rows ||
    !"parameter_rank" %in% names(overall)) {
  stop("The parameter and method summaries do not form a complete six-method grid.")
}
interval_columns <- c(
  "replicates", "ci95_lower_normalized_rf", "ci95_upper_normalized_rf"
)
if (any(!interval_columns %in% names(overall)) ||
    any(!interval_columns %in% names(methods)) ||
    any(!interval_columns %in% names(by_shape)) ||
    any(!interval_columns %in% names(methods_by_shape))) {
  stop("Run the current summarize-recorder-parameter-grid.R before plotting.")
}
replicate_counts <- unique(overall$replicates)
if (length(replicate_counts) != 1L || replicate_counts < 2L) {
  stop("At least two population-seed replicates are required for error bars.")
}
replicate_count <- replicate_counts[[1L]]
if (nrow(by_shape) != nrow(overall) * 4L ||
    nrow(methods_by_shape) != nrow(by_shape) * 6L ||
    !"shape" %in% names(by_shape)) {
  stop("The shape summaries do not form a complete four-shape, six-method grid.")
}

recorder_levels <- c(
  "baseline", "wt_crispr", "prime", "palincode", "mitochondrial"
)
# Facet labels list the rates actually present in the summary, in tier order,
# so an extended run does not advertise the original five values.
recorder_titles <- c(
  baseline = "BASELINE",
  wt_crispr = "WT CRISPR (FLARE-like)",
  prime = "Prime editing",
  palincode = "PALINCODE",
  mitochondrial = "Mitochondrial"
)
recorder_units <- c(
  baseline = "", wt_crispr = " ins/del", prime = "",
  palincode = "", mitochondrial = ""
)
format_tested_values <- function(system_name) {
  rows <- overall[overall$system == system_name, , drop = FALSE]
  if (!nrow(rows)) return("")
  rows <- rows[order(rows$parameter_rank), , drop = FALSE]
  primary <- formatC(rows$parameter_value, format = "g", digits = 4)
  formatted <- if (all(is.na(rows$parameter_value_secondary))) {
    primary
  } else {
    sprintf(
      "(%s, %s)", primary,
      formatC(rows$parameter_value_secondary, format = "g", digits = 4)
    )
  }
  paste(formatted, collapse = " / ")
}
recorder_labels <- vapply(recorder_levels, function(system_name) {
  paste0(
    recorder_titles[[system_name]], "\n",
    format_tested_values(system_name), recorder_units[[system_name]]
  )
}, character(1))
names(recorder_labels) <- recorder_levels
recorder_short_labels <- c(
  baseline = "BASELINE",
  wt_crispr = "WT CRISPR\n(FLARE-like)",
  prime = "Prime editing",
  palincode = "PALINCODE",
  mitochondrial = "Mitochondrial"
)
# Tier names and their axis order come from the summaries rather than a fixed
# five-tier list, so the extension tiers above very_high (very_high_2x and
# beyond, added to bracket the BASELINE and mitochondrial optima) plot without
# further changes here.
observed_settings <- unique(rbind(
  overall[c("parameter_setting", "parameter_rank")],
  methods[c("parameter_setting", "parameter_rank")],
  by_shape[c("parameter_setting", "parameter_rank")],
  methods_by_shape[c("parameter_setting", "parameter_rank")]
))
observed_settings <- observed_settings[
  order(observed_settings$parameter_rank), , drop = FALSE
]
setting_levels <- as.character(observed_settings$parameter_setting)
setting_breaks <- observed_settings$parameter_rank
setting_limits <- c(min(setting_breaks) - 0.35, max(setting_breaks) + 0.35)

# "very_high_4x" is four times the very_high rate; the base tiers keep their
# original wording.
setting_long_labels <- function(value) {
  base <- c(
    very_low = "Very low", low = "Low", mid = "Mid",
    high = "High", very_high = "Very high"
  )
  labels <- unname(base[value])
  multiple <- sub("^very_high_([0-9]+)x$", "\\1", value)
  extended <- grepl("^very_high_[0-9]+x$", value)
  labels[extended] <- paste0("Very high\n×", multiple[extended])
  labels
}
setting_short_labels <- function(value) {
  base <- c(
    very_low = "VL", low = "L", mid = "M", high = "H", very_high = "VH"
  )
  labels <- unname(base[value])
  multiple <- sub("^very_high_([0-9]+)x$", "\\1", value)
  extended <- grepl("^very_high_[0-9]+x$", value)
  labels[extended] <- paste0("VH×", multiple[extended])
  labels
}
shape_levels <- c("balanced", "comb", "hierarchical", "neutral")
shape_labels <- c(
  balanced = "Balanced",
  comb = "Comb",
  hierarchical = "Hierarchical",
  neutral = "Neutral"
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

prepare_data <- function(value) {
  value$system <- factor(value$system, levels = recorder_levels)
  value$parameter_setting <- factor(
    value$parameter_setting, levels = setting_levels
  )
  value$setting_index <- value$parameter_rank
  value
}
overall <- prepare_data(overall)
methods <- prepare_data(methods)
by_shape <- prepare_data(by_shape)
methods_by_shape <- prepare_data(methods_by_shape)
methods$method <- factor(methods$method, levels = method_levels)
by_shape$shape <- factor(by_shape$shape, levels = shape_levels)
methods_by_shape$shape <- factor(
  methods_by_shape$shape, levels = shape_levels
)
methods_by_shape$method <- factor(
  methods_by_shape$method, levels = method_levels
)
method_offsets <- stats::setNames(
  seq(-0.10, 0.10, length.out = length(method_levels)), method_levels
)
methods$plot_x <- methods$setting_index + method_offsets[as.character(methods$method)]
methods_by_shape$plot_x <- methods_by_shape$setting_index +
  method_offsets[as.character(methods_by_shape$method)]
if (anyNA(by_shape$shape) || anyNA(methods_by_shape$shape)) {
  stop("The shape summaries contain an unexpected tree shape.")
}

best_rows <- do.call(rbind, lapply(
  split(overall, overall$system),
  function(part) part[which.min(part$mean_normalized_rf), , drop = FALSE]
))
rownames(best_rows) <- NULL

shape_best_rows <- do.call(rbind, lapply(
  split(by_shape, interaction(by_shape$system, by_shape$shape, drop = TRUE)),
  function(part) part[which.min(part$mean_normalized_rf), , drop = FALSE]
))
rownames(shape_best_rows) <- NULL

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 12, base_family = "sans"))
plot_value <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = best_rows,
    ggplot2::aes(
      xmin = setting_index - 0.34,
      xmax = setting_index + 0.34,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "#E9C75C",
    alpha = 0.20,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_errorbar(
    data = methods,
    ggplot2::aes(
      x = plot_x,
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf,
      color = method
    ),
    width = 0.055,
    linewidth = 0.55,
    alpha = 0.68
  ) +
  ggplot2::geom_line(
    data = methods,
    ggplot2::aes(
      x = plot_x,
      y = mean_normalized_rf,
      color = method,
      group = method
    ),
    linewidth = 0.75,
    alpha = 0.78
  ) +
  ggplot2::geom_point(
    data = methods,
    ggplot2::aes(
      x = plot_x,
      y = mean_normalized_rf,
      color = method
    ),
    size = 2.25,
    alpha = 0.92
  ) +
  ggplot2::geom_errorbar(
    data = overall,
    ggplot2::aes(
      x = setting_index,
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf
    ),
    width = 0.15,
    linewidth = 0.8,
    color = "#171717"
  ) +
  ggplot2::geom_line(
    data = overall,
    ggplot2::aes(
      x = setting_index,
      y = mean_normalized_rf,
      group = system
    ),
    color = "#171717",
    linewidth = 1.25
  ) +
  ggplot2::geom_point(
    data = overall,
    ggplot2::aes(x = setting_index, y = mean_normalized_rf),
    shape = 21,
    size = 3.5,
    stroke = 1.0,
    color = "#171717",
    fill = "#FFFDF7"
  ) +
  ggplot2::geom_point(
    data = best_rows,
    ggplot2::aes(x = setting_index, y = mean_normalized_rf),
    shape = 23,
    size = 4.6,
    stroke = 1.0,
    color = "#171717",
    fill = "#E9C75C"
  ) +
  ggplot2::geom_text(
    data = overall,
    ggplot2::aes(
      x = setting_index,
      y = mean_normalized_rf,
      label = sprintf("%.3f", mean_normalized_rf)
    ),
    color = "#171717",
    size = 3.1,
    fontface = "bold",
    vjust = -1.15
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(system),
    ncol = 2,
    labeller = ggplot2::as_labeller(recorder_labels)
  ) +
  ggplot2::scale_x_continuous(
    breaks = setting_breaks,
    labels = setting_long_labels(setting_levels),
    limits = setting_limits
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0.3, 1.0, by = 0.1),
    limits = c(0.30, 1.02),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_color_manual(
    values = method_colors,
    breaks = method_levels,
    labels = method_labels
  ) +
  ggplot2::labs(
    title = "Editing-rate sensitivity across lineage recorders",
    subtitle = paste(
      "Normalized Robinson-Foulds distance across four tree shapes;",
      "lower values indicate better reconstruction"
    ),
    x = "Recorder-specific editing-rate setting",
    y = "Mean normalized RF distance",
    color = "Tree-building method",
    caption = paste(
      "Colored lines: method means. Black line: mean across all six methods.",
      "Error bars: 95% CIs across", replicate_count,
      "population seeds. Gold marks: best tested setting. n = 250 cells."
    )
  ) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "#FFFFFF", color = NA),
    panel.background = ggplot2::element_rect(fill = "#FFFFFF", color = NA),
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      fill = "#FFFFFF", color = "#4F5D5F", linewidth = 0.7
    ),
    strip.text = ggplot2::element_text(
      color = "#384849", face = "bold", size = 11, lineheight = 1.1,
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
    legend.key.width = grid::unit(1.35, "lines"),
    panel.spacing = grid::unit(1.0, "lines"),
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
  height = 11.2,
  units = "in",
  dpi = 300,
  background = "#FFFFFF"
)
ggplot2::ggsave(
  pdf_path,
  plot = plot_value,
  device = grDevices::cairo_pdf,
  width = 12,
  height = 11.2,
  units = "in",
  bg = "#FFFFFF"
)

shape_plot <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = shape_best_rows,
    ggplot2::aes(
      xmin = setting_index - 0.34,
      xmax = setting_index + 0.34,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "#E9C75C",
    alpha = 0.20,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_errorbar(
    data = methods_by_shape,
    ggplot2::aes(
      x = plot_x,
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf,
      color = method
    ),
    width = 0.045,
    linewidth = 0.45,
    alpha = 0.65
  ) +
  ggplot2::geom_line(
    data = methods_by_shape,
    ggplot2::aes(
      x = plot_x,
      y = mean_normalized_rf,
      color = method,
      group = method
    ),
    linewidth = 0.62,
    alpha = 0.78
  ) +
  ggplot2::geom_point(
    data = methods_by_shape,
    ggplot2::aes(
      x = plot_x,
      y = mean_normalized_rf,
      color = method
    ),
    size = 1.45,
    alpha = 0.90
  ) +
  ggplot2::geom_errorbar(
    data = by_shape,
    ggplot2::aes(
      x = setting_index,
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf
    ),
    width = 0.13,
    linewidth = 0.72,
    color = "#171717"
  ) +
  ggplot2::geom_line(
    data = by_shape,
    ggplot2::aes(
      x = setting_index,
      y = mean_normalized_rf,
      group = interaction(system, shape)
    ),
    color = "#171717",
    linewidth = 1.05
  ) +
  ggplot2::geom_point(
    data = by_shape,
    ggplot2::aes(x = setting_index, y = mean_normalized_rf),
    shape = 21,
    size = 2.5,
    stroke = 0.85,
    color = "#171717",
    fill = "#FFFDF7"
  ) +
  ggplot2::geom_point(
    data = shape_best_rows,
    ggplot2::aes(x = setting_index, y = mean_normalized_rf),
    shape = 23,
    size = 3.5,
    stroke = 0.9,
    color = "#171717",
    fill = "#E9C75C"
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(shape),
    cols = ggplot2::vars(system),
    labeller = ggplot2::labeller(
      shape = ggplot2::as_labeller(shape_labels),
      system = ggplot2::as_labeller(recorder_short_labels)
    ),
    switch = "y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = setting_breaks,
    labels = setting_short_labels(setting_levels),
    limits = setting_limits
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1.02),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_color_manual(
    values = method_colors,
    breaks = method_levels,
    labels = method_labels
  ) +
  ggplot2::labs(
    title = "Editing-rate sensitivity by tree shape",
    subtitle = paste(
      "Each panel separates one recorder and population topology;",
      "lower normalized Robinson-Foulds distance is better"
    ),
    x = "Recorder-specific setting (very low to very high)",
    y = "Mean normalized RF distance",
    color = "Tree-building method",
    caption = paste(
      "Colored lines: individual methods. Black line: mean across all six methods.",
      "Error bars: 95% CIs across", replicate_count,
      "population seeds. Gold marks: best setting within each panel. n = 250."
    )
  ) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "#FFFFFF", color = NA),
    panel.background = ggplot2::element_rect(fill = "#FFFFFF", color = NA),
    panel.grid = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(
      fill = "#FFFFFF", color = "#4F5D5F", linewidth = 0.7
    ),
    strip.text = ggplot2::element_text(
      color = "#384849", face = "bold", size = 10,
      margin = ggplot2::margin(7, 6, 7, 6)
    ),
    strip.text.y.left = ggplot2::element_text(angle = 0),
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
    legend.key.width = grid::unit(1.25, "lines"),
    panel.spacing = grid::unit(0.65, "lines"),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

shape_png_path <- paste0(shape_output_prefix, ".png")
shape_pdf_path <- paste0(shape_output_prefix, ".pdf")
dir.create(dirname(shape_png_path), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  shape_png_path,
  plot = shape_plot,
  device = ragg::agg_png,
  width = 18,
  height = 13,
  units = "in",
  dpi = 300,
  background = "#FFFFFF"
)
ggplot2::ggsave(
  shape_pdf_path,
  plot = shape_plot,
  device = grDevices::cairo_pdf,
  width = 18,
  height = 13,
  units = "in",
  bg = "#FFFFFF"
)

neutral_overall <- by_shape[as.character(by_shape$shape) == "neutral", , drop = FALSE]
neutral_methods <- methods_by_shape[
  as.character(methods_by_shape$shape) == "neutral", , drop = FALSE
]
neutral_best_rows <- shape_best_rows[
  as.character(shape_best_rows$shape) == "neutral", , drop = FALSE
]
if (nrow(neutral_overall) != nrow(overall) ||
    nrow(neutral_methods) != nrow(methods) ||
    nrow(neutral_best_rows) != length(recorder_levels)) {
  stop("The neutral-tree summaries do not form a complete recorder grid.")
}

neutral_plot <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = neutral_best_rows,
    ggplot2::aes(
      xmin = setting_index - 0.34,
      xmax = setting_index + 0.34,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "#E9C75C",
    alpha = 0.20,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_errorbar(
    data = neutral_methods,
    ggplot2::aes(
      x = plot_x,
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf,
      color = method
    ),
    width = 0.045,
    linewidth = 0.45,
    alpha = 0.65
  ) +
  ggplot2::geom_line(
    data = neutral_methods,
    ggplot2::aes(
      x = plot_x,
      y = mean_normalized_rf,
      color = method,
      group = method
    ),
    linewidth = 0.68,
    alpha = 0.80
  ) +
  ggplot2::geom_point(
    data = neutral_methods,
    ggplot2::aes(
      x = plot_x,
      y = mean_normalized_rf,
      color = method
    ),
    size = 1.8,
    alpha = 0.92
  ) +
  ggplot2::geom_errorbar(
    data = neutral_overall,
    ggplot2::aes(
      x = setting_index,
      ymin = ci95_lower_normalized_rf,
      ymax = ci95_upper_normalized_rf
    ),
    width = 0.13,
    linewidth = 0.76,
    color = "#171717"
  ) +
  ggplot2::geom_line(
    data = neutral_overall,
    ggplot2::aes(
      x = setting_index,
      y = mean_normalized_rf,
      group = system
    ),
    color = "#171717",
    linewidth = 1.12
  ) +
  ggplot2::geom_point(
    data = neutral_overall,
    ggplot2::aes(x = setting_index, y = mean_normalized_rf),
    shape = 21,
    size = 2.8,
    stroke = 0.9,
    color = "#171717",
    fill = "#FFFDF7"
  ) +
  ggplot2::geom_point(
    data = neutral_best_rows,
    ggplot2::aes(x = setting_index, y = mean_normalized_rf),
    shape = 23,
    size = 3.8,
    stroke = 0.95,
    color = "#171717",
    fill = "#E9C75C"
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(system),
    nrow = 1,
    labeller = ggplot2::as_labeller(recorder_short_labels)
  ) +
  ggplot2::scale_x_continuous(
    breaks = setting_breaks,
    labels = setting_short_labels(setting_levels),
    limits = setting_limits
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1.02),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_color_manual(
    values = method_colors,
    breaks = method_levels,
    labels = method_labels
  ) +
  ggplot2::labs(
    title = "Editing-rate sensitivity for neutral lineage trees",
    subtitle = paste(
      "Pure-birth Gillespie populations with 250 sampled cells;",
      "lower normalized Robinson-Foulds distance is better"
    ),
    x = "Recorder-specific setting (very low to very high)",
    y = "Mean normalized RF distance",
    color = "Tree-building method",
    caption = paste(
      "Colored lines: individual methods. Black line: mean across all six methods.",
      "Error bars: 95% CIs across", replicate_count,
      "population seeds. Gold marks: best neutral-tree setting."
    )
  ) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "#FFFFFF", color = NA),
    panel.background = ggplot2::element_rect(fill = "#FFFFFF", color = NA),
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      fill = "#FFFFFF", color = "#4F5D5F", linewidth = 0.7
    ),
    strip.text = ggplot2::element_text(
      color = "#384849", face = "bold", size = 10,
      margin = ggplot2::margin(7, 6, 7, 6)
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
    legend.key.width = grid::unit(1.25, "lines"),
    panel.spacing = grid::unit(0.70, "lines"),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

neutral_png_path <- paste0(neutral_output_prefix, ".png")
neutral_pdf_path <- paste0(neutral_output_prefix, ".pdf")
dir.create(dirname(neutral_png_path), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  neutral_png_path,
  plot = neutral_plot,
  device = ragg::agg_png,
  width = 18,
  height = 5.8,
  units = "in",
  dpi = 300,
  background = "#FFFFFF"
)
ggplot2::ggsave(
  neutral_pdf_path,
  plot = neutral_plot,
  device = grDevices::cairo_pdf,
  width = 18,
  height = 5.8,
  units = "in",
  bg = "#FFFFFF"
)
cat("PNG: ", normalizePath(png_path), "\n", sep = "")
cat("PDF: ", normalizePath(pdf_path), "\n", sep = "")
cat("Shape PNG: ", normalizePath(shape_png_path), "\n", sep = "")
cat("Shape PDF: ", normalizePath(shape_pdf_path), "\n", sep = "")
cat("Neutral PNG: ", normalizePath(neutral_png_path), "\n", sep = "")
cat("Neutral PDF: ", normalizePath(neutral_pdf_path), "\n", sep = "")
