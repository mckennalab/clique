#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1]] else {
  "test_harness/tree_building/results/prime_rate_grid_n250_k05"
}
root <- normalizePath(root, mustWork = TRUE)
simulation_dir <- file.path(root, "simulation")
reconstruction_dir <- file.path(root, "reconstruction")

manifest <- read.csv(
  gzfile(file.path(simulation_dir, "benchmark_manifest.csv.gz")),
  check.names = FALSE
)
results <- read.csv(
  gzfile(file.path(reconstruction_dir, "reconstruction_summary.csv.gz")),
  check.names = FALSE
)
results$editing_rate <- manifest$editing_rate[
  match(results$condition_id, manifest$condition_id)
]
if (anyNA(results$editing_rate)) {
  stop("A reconstruction condition is absent from the rate-grid manifest.")
}

summarize_values <- function(data, groups) {
  keys <- interaction(data[groups], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(data, keys), function(part) {
    values <- part$normalized_rf[part$status == "success"]
    row <- part[1L, groups, drop = FALSE]
    row$conditions <- nrow(part)
    row$successful <- length(values)
    row$mean_normalized_rf <- if (length(values)) mean(values) else NA_real_
    row$median_normalized_rf <- if (length(values)) median(values) else NA_real_
    row
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output[do.call(order, unname(output[groups])), , drop = FALSE]
}

rate_method <- summarize_values(results, c("editing_rate", "method"))
rate_shape_method <- summarize_values(
  results, c("editing_rate", "shape", "method")
)
rate_overall <- summarize_values(results, "editing_rate")
rate_shape <- summarize_values(results, c("editing_rate", "shape"))
best_rate <- do.call(rbind, lapply(split(rate_method, rate_method$method), function(part) {
  part[which.min(part$mean_normalized_rf), c(
    "method", "editing_rate", "mean_normalized_rf", "median_normalized_rf"
  )]
}))
rownames(best_rate) <- NULL
best_rate <- best_rate[order(best_rate$mean_normalized_rf), , drop = FALSE]

split_key <- function(labels, all_labels) {
  side_a <- sort(labels)
  side_b <- sort(setdiff(all_labels, side_a))
  key_a <- paste(side_a, collapse = "|")
  key_b <- paste(side_b, collapse = "|")
  if (length(side_a) < length(side_b) ||
      (length(side_a) == length(side_b) && key_a <= key_b)) key_a else key_b
}

tree_splits <- function(tree) {
  number_tips <- length(tree$tip.label)
  partitions <- ape::prop.part(tree)
  retained <- lengths(partitions) >= 2L & lengths(partitions) <= number_tips - 2L
  unique(vapply(partitions[retained], function(descendants) {
    split_key(tree$tip.label[descendants], tree$tip.label)
  }, character(1)))
}

truth_splits_by_shape <- new.env(parent = emptyenv())
diagnostics <- lapply(seq_len(nrow(manifest)), function(index) {
  matrix_value <- as.matrix(readRDS(file.path(
    simulation_dir, manifest$character_matrix[index]
  )))
  shape <- manifest$shape[index]
  if (!exists(shape, envir = truth_splits_by_shape, inherits = FALSE)) {
    truth <- ape::read.tree(file.path(
      simulation_dir, manifest$ground_truth_tree[index]
    ))
    assign(shape, tree_splits(truth), envir = truth_splits_by_shape)
  }
  prevalence <- colMeans(matrix_value != 0)
  edited_counts <- colSums(matrix_value != 0)
  informative <- which(
    edited_counts >= 2L & edited_counts <= nrow(matrix_value) - 2L
  )
  truth_keys <- get(shape, envir = truth_splits_by_shape, inherits = FALSE)
  exact <- vapply(informative, function(column) {
    split_key(
      rownames(matrix_value)[matrix_value[, column] != 0],
      rownames(matrix_value)
    ) %in% truth_keys
  }, logical(1))
  data.frame(
    shape = shape,
    editing_rate = manifest$editing_rate[index],
    cells = nrow(matrix_value),
    characters = ncol(matrix_value),
    unique_profiles = nrow(unique(matrix_value)),
    variable_characters = sum(prevalence > 0 & prevalence < 1),
    informative_characters = length(informative),
    exact_true_split_characters = sum(exact),
    edited_fraction = mean(matrix_value != 0),
    median_site_prevalence = median(prevalence),
    max_site_prevalence = max(prevalence)
  )
})
diagnostics <- do.call(rbind, diagnostics)

write_tsv <- function(value, name) {
  write.table(
    value, file.path(root, name), sep = "\t", row.names = FALSE,
    quote = FALSE
  )
}
write_tsv(rate_method, "accuracy_by_rate_method.tsv")
write_tsv(rate_shape_method, "accuracy_by_rate_shape_method.tsv")
write_tsv(rate_overall, "accuracy_by_rate_overall.tsv")
write_tsv(rate_shape, "accuracy_by_rate_shape.tsv")
write_tsv(best_rate, "best_rate_by_method.tsv")
write_tsv(diagnostics, "matrix_diagnostics.tsv")

cat("Accuracy by editing rate:\n")
print(rate_overall, row.names = FALSE)
cat("\nBest editing rate by method:\n")
print(best_rate, row.names = FALSE)
