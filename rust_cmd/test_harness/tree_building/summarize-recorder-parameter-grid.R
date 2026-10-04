#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else {
  "test_harness/tree_building/results/recorder_parameter_grid_n250_level05"
}
root <- normalizePath(root, mustWork = TRUE)

simulation_dir <- file.path(root, "simulation")
reconstruction_path <- file.path(
  root, "reconstruction", "reconstruction_summary.csv.gz"
)
manifest_path <- file.path(simulation_dir, "benchmark_manifest.csv.gz")
diagnostics_path <- file.path(simulation_dir, "matrix_diagnostics.tsv")
required <- c(reconstruction_path, manifest_path, diagnostics_path)
if (any(!file.exists(required))) {
  stop("Missing parameter-grid input: ", paste(required[!file.exists(required)], collapse = ", "))
}

results <- read.csv(gzfile(reconstruction_path), check.names = FALSE)
manifest <- read.csv(gzfile(manifest_path), check.names = FALSE)
diagnostics <- read.delim(diagnostics_path, check.names = FALSE)
# Reconstruction methods the run is required to contain, as a comma-separated
# second argument. The default is the full six-method panel used by the
# published grid; the bracket search passes a smaller set (typically just nj)
# because MIX and IQ-TREE cost hours per tier and are not needed to locate a
# minimum. The completeness check below is enforced against whichever set is
# named, so a partial run still cannot be summarized by accident.
methods <- if (length(arguments) >= 2L) {
  requested <- trimws(strsplit(arguments[[2L]], ",", fixed = TRUE)[[1L]])
  requested[nzchar(requested)]
} else {
  c("nj", "parsimony", "iqtree", "mix", "vine", "cassiopeia")
}
if (!length(methods) || anyDuplicated(methods)) {
  stop("The reconstruction method list must be non-empty and unique.")
}

if (!nrow(manifest) || length(unique(manifest$condition_id)) != nrow(manifest) ||
    !"parameter_rank" %in% names(manifest)) {
  stop("The simulation manifest must contain unique ranked parameter conditions.")
}
if (nrow(results) != nrow(manifest) * length(methods) ||
    any(results$status != "success") ||
    !setequal(unique(results$method), methods) ||
    anyDuplicated(results[c("condition_id", "method")])) {
  stop("The reconstruction grid is incomplete or contains failed/duplicate results.")
}

parameter_columns <- c(
  "parameter_setting", "parameter_rank", "editing_parameter", "parameter_value",
  "parameter_value_secondary"
)
manifest_index <- match(results$condition_id, manifest$condition_id)
if (anyNA(manifest_index)) {
  stop("A reconstructed condition is absent from the simulation manifest.")
}
results <- cbind(
  results,
  manifest[manifest_index, parameter_columns, drop = FALSE]
)

summarize_groups <- function(data, groups) {
  key_data <- lapply(data[groups], function(value) {
    value <- as.character(value)
    value[is.na(value)] <- "<NA>"
    value
  })
  keys <- interaction(key_data, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(data, keys), function(part) {
    replicate_means <- vapply(
      split(part$normalized_rf, part$population_seed),
      mean,
      numeric(1)
    )
    replicate_count <- length(replicate_means)
    replicate_sd <- if (replicate_count > 1L) {
      stats::sd(replicate_means)
    } else {
      NA_real_
    }
    replicate_se <- replicate_sd / sqrt(replicate_count)
    confidence_radius <- if (replicate_count > 1L) {
      stats::qt(0.975, df = replicate_count - 1L) * replicate_se
    } else {
      NA_real_
    }
    row <- part[1L, groups, drop = FALSE]
    row$trees <- nrow(part)
    row$replicates <- replicate_count
    row$trees_per_replicate_min <- min(table(part$population_seed))
    row$trees_per_replicate_max <- max(table(part$population_seed))
    row$mean_normalized_rf <- mean(replicate_means)
    row$replicate_sd_normalized_rf <- replicate_sd
    row$replicate_se_normalized_rf <- replicate_se
    row$ci95_lower_normalized_rf <- max(
      0, row$mean_normalized_rf - confidence_radius, na.rm = TRUE
    )
    row$ci95_upper_normalized_rf <- min(
      1, row$mean_normalized_rf + confidence_radius, na.rm = TRUE
    )
    if (is.na(confidence_radius)) {
      row$ci95_lower_normalized_rf <- NA_real_
      row$ci95_upper_normalized_rf <- NA_real_
    }
    row$median_normalized_rf <- median(part$normalized_rf)
    row$minimum_normalized_rf <- min(part$normalized_rf)
    row$maximum_normalized_rf <- max(part$normalized_rf)
    row$mean_elapsed_seconds <- mean(part$elapsed_seconds)
    row
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output[do.call(order, unname(output[groups])), , drop = FALSE]
}

parameter_groups <- c(
  "system", "parameter_setting", "parameter_rank", "editing_parameter", "parameter_value",
  "parameter_value_secondary"
)
by_parameter <- summarize_groups(results, parameter_groups)
by_parameter_method <- summarize_groups(
  results, c(parameter_groups, "method")
)
by_parameter_shape <- summarize_groups(
  results, c(parameter_groups, "shape")
)
by_parameter_shape_method <- summarize_groups(
  results, c(parameter_groups, "shape", "method")
)

best_parameter <- do.call(rbind, lapply(
  split(by_parameter, by_parameter$system),
  function(part) part[which.min(part$mean_normalized_rf), , drop = FALSE]
))
rownames(best_parameter) <- NULL
best_parameter <- best_parameter[
  order(best_parameter$mean_normalized_rf), , drop = FALSE
]

best_method <- do.call(rbind, lapply(
  split(by_parameter_method, interaction(
    by_parameter_method$system,
    by_parameter_method$parameter_setting,
    drop = TRUE
  )),
  function(part) part[which.min(part$mean_normalized_rf), c(
    parameter_groups, "method", "mean_normalized_rf", "median_normalized_rf"
  ), drop = FALSE]
))
rownames(best_method) <- NULL
best_method <- best_method[
  do.call(order, unname(best_method[c("system", "parameter_rank")])),
  , drop = FALSE
]

best_parameter_method <- do.call(rbind, lapply(
  split(by_parameter_method, interaction(
    by_parameter_method$system,
    by_parameter_method$method,
    drop = TRUE
  )),
  function(part) part[which.min(part$mean_normalized_rf), c(
    parameter_groups, "method", "mean_normalized_rf", "median_normalized_rf",
    "mean_elapsed_seconds"
  ), drop = FALSE]
))
rownames(best_parameter_method) <- NULL
best_parameter_method <- best_parameter_method[
  do.call(order, unname(best_parameter_method[c("system", "method")])),
  , drop = FALSE
]

diagnostic_groups <- c(
  "system", "parameter_setting", "parameter_rank", "parameter_value",
  "parameter_value_secondary"
)
diagnostic_keys <- interaction(
  lapply(diagnostics[diagnostic_groups], function(value) {
    value <- as.character(value)
    value[is.na(value)] <- "<NA>"
    value
  }),
  drop = TRUE, lex.order = TRUE
)
diagnostic_summary <- do.call(rbind, lapply(
  split(diagnostics, diagnostic_keys),
  function(part) {
    row <- part[1L, diagnostic_groups, drop = FALSE]
    row$conditions <- nrow(part)
    row$variable_characters_min <- min(part$variable_characters)
    row$variable_characters_max <- max(part$variable_characters)
    row$informative_characters_min <- min(part$informative_characters)
    row$informative_characters_max <- max(part$informative_characters)
    row$unique_profiles_min <- min(part$unique_profiles)
    row$unique_profiles_max <- max(part$unique_profiles)
    row$mean_nonzero_fraction <- mean(part$nonzero_fraction)
    row
  }
))
rownames(diagnostic_summary) <- NULL
diagnostic_summary <- diagnostic_summary[
  do.call(order, unname(diagnostic_summary[c("system", "parameter_rank")])),
  , drop = FALSE
]

write_tsv <- function(value, name) {
  write.table(
    value, file.path(root, name), sep = "\t", row.names = FALSE,
    quote = FALSE
  )
}

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The ape package is required to validate reconstructed trees.")
}
resolve_output_path <- function(path, base_dir) {
  if (file.exists(path)) path else file.path(base_dir, path)
}
truth_tip_sets <- lapply(seq_len(nrow(manifest)), function(index) {
  truth_path <- resolve_output_path(
    manifest$ground_truth_tree[[index]], simulation_dir
  )
  sort(ape::read.tree(truth_path)$tip.label)
})
names(truth_tip_sets) <- manifest$condition_id
validation_rows <- lapply(seq_len(nrow(results)), function(index) {
  tree_path <- resolve_output_path(
    results$output_tree[[index]], file.path(root, "reconstruction")
  )
  expected_tips <- truth_tip_sets[[results$condition_id[[index]]]]
  tryCatch({
    tree <- ape::read.tree(tree_path)
    observed_tips <- sort(tree$tip.label)
    data.frame(
      condition_id = results$condition_id[[index]],
      method = results$method[[index]],
      tips = length(observed_tips),
      expected_tips = length(expected_tips),
      readable = TRUE,
      exact_tip_set = identical(observed_tips, expected_tips),
      error = ""
    )
  }, error = function(error) {
    data.frame(
      condition_id = results$condition_id[[index]],
      method = results$method[[index]],
      tips = NA_integer_,
      expected_tips = length(expected_tips),
      readable = FALSE,
      exact_tip_set = FALSE,
      error = conditionMessage(error)
    )
  })
})
tree_validation <- do.call(rbind, validation_rows)
if (any(!tree_validation$readable) || any(!tree_validation$exact_tip_set)) {
  stop("At least one reconstructed tree failed Newick/tip-set validation.")
}

write_tsv(by_parameter, "accuracy_by_recorder_parameter.tsv")
write_tsv(by_parameter_method, "accuracy_by_recorder_parameter_method.tsv")
write_tsv(by_parameter_shape, "accuracy_by_recorder_parameter_shape.tsv")
write_tsv(
  by_parameter_shape_method,
  "accuracy_by_recorder_parameter_shape_method.tsv"
)
write_tsv(best_parameter, "best_parameter_by_recorder.tsv")
write_tsv(best_method, "best_method_by_recorder_parameter.tsv")
write_tsv(
  best_parameter_method,
  "best_parameter_by_recorder_method.tsv"
)
write_tsv(diagnostic_summary, "matrix_diagnostics_by_recorder_parameter.tsv")
write_tsv(tree_validation, "tree_validation.tsv")

joined_path <- gzfile(
  file.path(root, "reconstruction_with_parameters.tsv.gz"), open = "wt"
)
write.table(
  results, joined_path, sep = "\t", row.names = FALSE, quote = FALSE
)
close(joined_path)

cat("Accuracy by recorder parameter:\n")
print(
  by_parameter[order(by_parameter$system, by_parameter$parameter_rank), c(
    "system", "parameter_setting", "parameter_value",
    "parameter_value_secondary", "mean_normalized_rf", "median_normalized_rf"
  )],
  row.names = FALSE
)
cat("\nBest parameter by recorder:\n")
print(best_parameter, row.names = FALSE)
