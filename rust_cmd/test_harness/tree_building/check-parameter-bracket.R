#!/usr/bin/env Rscript

# Reports whether each recorder's editing-rate optimum is bracketed by the
# tested grid.
#
# An optimum is bracketed when the tier with the lowest mean normalized RF has
# at least one tested tier on each side of it. When the best tier sits at either
# end of the ladder, the true optimum may lie beyond the grid and the reported
# best setting is only "the best of what was tried" — the situation this check
# exists to detect. BASELINE and mitochondrial recording both selected the top
# tier in the original five-tier grid, which is why the grid is extended upward.
#
# Bracketing is reported both across all methods present (from
# accuracy_by_recorder_parameter.tsv) and per reconstruction method (from
# accuracy_by_recorder_parameter_method.tsv), because the two can disagree: a
# method-specific optimum can sit on the edge while the across-method mean is
# interior.
#
# A bracket is also called weak when the neighbouring tier's 95% confidence
# interval contains the best tier's mean, i.e. the minimum is not separated
# from its neighbour by more than replicate noise. A weak bracket is a real
# bracket in shape but should not be quoted as a resolved optimum.
#
# Usage:
#   Rscript check-parameter-bracket.R <results-root> [output.tsv]
#
# Writes bracket_status.tsv under the results root (or the named path) and
# prints the same table. Exit status is 0 when every recorder is bracketed
# across methods and 1 when at least one is not, so a shell loop can use it as
# a stopping condition.

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else {
  "test_harness/tree_building/results/recorder_parameter_grid_n250_level05"
}
root <- normalizePath(root, mustWork = TRUE)
output_path <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  file.path(root, "bracket_status.tsv")
}

parameter_path <- file.path(root, "accuracy_by_recorder_parameter.tsv")
method_path <- file.path(root, "accuracy_by_recorder_parameter_method.tsv")
if (any(!file.exists(c(parameter_path, method_path)))) {
  stop("Run summarize-recorder-parameter-grid.R before checking the bracket.")
}

overall <- read.delim(parameter_path, check.names = FALSE)
methods <- read.delim(method_path, check.names = FALSE)
required_columns <- c(
  "system", "parameter_setting", "parameter_rank", "parameter_value",
  "mean_normalized_rf"
)
if (any(!required_columns %in% names(overall)) ||
    any(!required_columns %in% names(methods))) {
  stop("The accuracy summaries are missing required columns.")
}

# Confidence bounds are optional so this also runs on a single-seed pilot,
# where the summarizer leaves them NA.
interval_columns <- c(
  "ci95_lower_normalized_rf", "ci95_upper_normalized_rf"
)
has_intervals <- all(interval_columns %in% names(overall)) &&
  all(interval_columns %in% names(methods))

# Describes one recorder's ladder: where the minimum sits, whether tiers exist
# on both sides of it, and whether the neighbouring tier is separated from it
# by more than replicate noise.
bracket_row <- function(part, label_columns) {
  part <- part[order(part$parameter_rank), , drop = FALSE]
  best_index <- which.min(part$mean_normalized_rf)
  best <- part[best_index, , drop = FALSE]
  lower_exists <- best_index > 1L
  upper_exists <- best_index < nrow(part)
  edge <- if (!lower_exists && !upper_exists) {
    "single_tier"
  } else if (!lower_exists) {
    "lower_edge"
  } else if (!upper_exists) {
    "upper_edge"
  } else {
    "interior"
  }

  # The better of the two immediate neighbours: the tier that would have to be
  # ruled out for the minimum to be called resolved.
  neighbour_indices <- c(best_index - 1L, best_index + 1L)
  neighbour_indices <- neighbour_indices[
    neighbour_indices >= 1L & neighbour_indices <= nrow(part)
  ]
  separated <- NA
  neighbour_setting <- NA_character_
  if (length(neighbour_indices)) {
    neighbours <- part[neighbour_indices, , drop = FALSE]
    closest <- neighbours[which.min(neighbours$mean_normalized_rf), , drop = FALSE]
    neighbour_setting <- as.character(closest$parameter_setting)
    if (has_intervals && !is.na(closest$ci95_lower_normalized_rf)) {
      separated <- closest$ci95_lower_normalized_rf >
        best$mean_normalized_rf
    }
  }

  row <- best[, label_columns, drop = FALSE]
  row$tiers_tested <- nrow(part)
  row$best_setting <- as.character(best$parameter_setting)
  row$best_rank <- best$parameter_rank
  row$best_parameter_value <- best$parameter_value
  row$best_mean_normalized_rf <- best$mean_normalized_rf
  row$lowest_rank_tested <- min(part$parameter_rank)
  row$highest_rank_tested <- max(part$parameter_rank)
  row$position <- edge
  row$bracketed <- edge == "interior"
  row$closest_neighbour_setting <- neighbour_setting
  row$separated_from_neighbour <- separated
  row$bracket_quality <- if (edge != "interior") {
    "unbracketed"
  } else if (isTRUE(separated)) {
    "bracketed"
  } else if (is.na(separated)) {
    "bracketed_unknown_separation"
  } else {
    "weak_bracket"
  }
  row$recommendation <- if (edge == "upper_edge") {
    "extend the ladder upward (raise --extension-doublings)"
  } else if (edge == "lower_edge") {
    "extend the ladder downward"
  } else if (edge == "single_tier") {
    "test more than one tier"
  } else if (identical(row$bracket_quality, "weak_bracket")) {
    "bracketed, but add seeds to separate it from its neighbour"
  } else {
    "bracketed"
  }
  row
}

overall_status <- do.call(rbind, lapply(
  split(overall, overall$system),
  bracket_row,
  label_columns = "system"
))
overall_status$method <- "all_methods_mean"
rownames(overall_status) <- NULL

method_status <- do.call(rbind, lapply(
  split(methods, interaction(methods$system, methods$method, drop = TRUE)),
  bracket_row,
  label_columns = c("system", "method")
))
rownames(method_status) <- NULL

status <- rbind(
  overall_status[c("system", "method", setdiff(names(overall_status), c("system", "method")))],
  method_status[c("system", "method", setdiff(names(method_status), c("system", "method")))]
)
status <- status[
  do.call(order, unname(status[c("system", "method")])), , drop = FALSE
]
rownames(status) <- NULL

write.table(
  status, output_path, sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Bracket status by recorder (mean across all methods present):\n")
print(
  overall_status[c(
    "system", "tiers_tested", "best_setting", "best_parameter_value",
    "best_mean_normalized_rf", "position", "bracket_quality", "recommendation"
  )],
  row.names = FALSE
)

unbracketed <- overall_status[!overall_status$bracketed, , drop = FALSE]
if (nrow(unbracketed)) {
  cat("\nStill unbracketed: ", paste(unbracketed$system, collapse = ", "), "\n",
      sep = "")
  cat("Wrote ", output_path, "\n", sep = "")
  quit(status = 1L)
}

cat("\nEvery recorder's optimum is bracketed.\n")
cat("Wrote ", output_path, "\n", sep = "")
