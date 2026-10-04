#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: finalize-neutral-best-integration-series.R OUTPUT_ROOT")
}
root <- normalizePath(arguments[[1L]], mustWork = TRUE)
accuracy_path <- file.path(root, "accuracy_by_recorder_integration.tsv")
settings_path <- file.path(root, "simulation", "benchmark_settings.rds")
if (!file.exists(accuracy_path) || !file.exists(settings_path)) {
  stop("The integration series is missing its accuracy table or settings.")
}

accuracy <- read.delim(accuracy_path, check.names = FALSE)
settings <- readRDS(settings_path)
levels <- as.integer(settings$integration_or_recovery_levels)
systems <- as.character(settings$systems)
if (!nrow(accuracy) || !length(levels) || !length(systems)) {
  stop("The integration series contains no completed conditions.")
}

# These duplicate the rates that generate-recorder-parameter-grid.R selects
# under --grid=best, so they have to be updated alongside its best_settings.
# They are labels only: the run itself takes its rate from the grid, so a stale
# entry here mislabels a correct result rather than changing it.
parameter_labels <- c(
  baseline = "0.00952",
  wt_crispr = "0.005 insertion / 0.05 deletion",
  # PEtracer, 3 edit sites x 8 marks, at the grid's 'high' tier.
  prime = "0.18",
  palincode = "0.08 / 0.05",
  mitochondrial = "8e-6"
)

rows <- lapply(systems, function(system) {
  part <- accuracy[accuracy$system == system, , drop = FALSE]
  if (!nrow(part)) {
    return(data.frame(
      system = system,
      best_editing_parameter = unname(parameter_labels[[system]]),
      tested_levels = "",
      stopping_status = "not_tested",
      stopping_level = NA_integer_,
      mean_normalized_rf = NA_real_,
      best_method = NA_character_,
      best_method_mean_normalized_rf = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  part <- part[order(part$level), , drop = FALSE]
  reached <- which(part$threshold_reached)
  final_index <- if (length(reached)) reached[[1L]] else nrow(part)
  final <- part[final_index, , drop = FALSE]
  status <- if (length(reached)) {
    "threshold_reached"
  } else if (max(part$level) == max(levels)) {
    "maximum_tested"
  } else {
    "incomplete"
  }
  data.frame(
    system = system,
    best_editing_parameter = unname(parameter_labels[[system]]),
    tested_levels = paste(part$level, collapse = ","),
    stopping_status = status,
    stopping_level = final$level,
    mean_normalized_rf = final$mean_normalized_rf,
    best_method = final$best_method,
    best_method_mean_normalized_rf = final$best_method_mean_normalized_rf,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, rows)
rownames(summary) <- NULL
write.table(
  summary,
  file.path(root, "stopping_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
print(summary, row.names = FALSE)
