#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else {
  "test_harness/tree_building/results/independent_recorders_n250_level05"
}
root <- normalizePath(root, mustWork = TRUE)
results <- read.csv(gzfile(file.path(
  root, "reconstruction", "reconstruction_summary.csv.gz"
)))
if (!nrow(results) || any(results$status != "success")) {
  stop("The independent recorder reconstruction grid is incomplete.")
}

summarize_groups <- function(data, groups) {
  keys <- interaction(data[groups], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(data, keys), function(part) {
    row <- part[1L, groups, drop = FALSE]
    row$conditions <- nrow(part)
    row$mean_normalized_rf <- mean(part$normalized_rf)
    row$median_normalized_rf <- median(part$normalized_rf)
    row$mean_elapsed_seconds <- mean(part$elapsed_seconds)
    row
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output[do.call(order, unname(output[groups])), , drop = FALSE]
}

recorder <- summarize_groups(results, "system")
recorder_method <- summarize_groups(results, c("system", "method"))
recorder_shape <- summarize_groups(results, c("system", "shape"))
recorder_shape_method <- summarize_groups(
  results, c("system", "shape", "method")
)
best_method <- do.call(rbind, lapply(
  split(recorder_method, recorder_method$system),
  function(part) part[which.min(part$mean_normalized_rf), c(
    "system", "method", "mean_normalized_rf", "median_normalized_rf"
  )]
))
rownames(best_method) <- NULL
best_method <- best_method[order(best_method$mean_normalized_rf), , drop = FALSE]

write_tsv <- function(value, name) {
  write.table(
    value, file.path(root, name), sep = "\t", row.names = FALSE,
    quote = FALSE
  )
}
write_tsv(recorder, "accuracy_by_recorder.tsv")
write_tsv(recorder_method, "accuracy_by_recorder_method.tsv")
write_tsv(recorder_shape, "accuracy_by_recorder_shape.tsv")
write_tsv(recorder_shape_method, "accuracy_by_recorder_shape_method.tsv")
write_tsv(best_method, "best_method_by_recorder.tsv")

rate_grid_path <- file.path(
  dirname(root), "prime_rate_grid_n250_k05", "reconstruction",
  "reconstruction_summary.csv.gz"
)
if (file.exists(rate_grid_path)) {
  prime_grid <- read.csv(gzfile(rate_grid_path))
  prime_rate_040 <- prime_grid[
    grepl("_prime_rate_040_", prime_grid$condition_id, fixed = TRUE) &
      prime_grid$status == "success",
    names(results),
    drop = FALSE
  ]
  expected_prime_rows <- length(unique(results$shape)) *
    length(unique(results$method))
  if (nrow(prime_rate_040) != expected_prime_rows) {
    stop("The rate-0.04 prime reconstruction panel is incomplete.")
  }
  calibrated <- rbind(
    results[results$system != "prime", , drop = FALSE],
    prime_rate_040
  )
  calibrated_recorder <- summarize_groups(calibrated, "system")
  calibrated_method <- summarize_groups(calibrated, c("system", "method"))
  calibrated_shape <- summarize_groups(calibrated, c("system", "shape"))
  calibrated_best <- do.call(rbind, lapply(
    split(calibrated_method, calibrated_method$system),
    function(part) part[which.min(part$mean_normalized_rf), c(
      "system", "method", "mean_normalized_rf", "median_normalized_rf"
    )]
  ))
  rownames(calibrated_best) <- NULL
  calibrated_best <- calibrated_best[
    order(calibrated_best$mean_normalized_rf), , drop = FALSE
  ]
  write_tsv(
    calibrated_recorder,
    "accuracy_by_recorder_prime_rate040.tsv"
  )
  write_tsv(
    calibrated_method,
    "accuracy_by_recorder_method_prime_rate040.tsv"
  )
  write_tsv(
    calibrated_shape,
    "accuracy_by_recorder_shape_prime_rate040.tsv"
  )
  write_tsv(
    calibrated_best,
    "best_method_by_recorder_prime_rate040.tsv"
  )
}

cat("Accuracy by independent recorder:\n")
print(recorder[order(recorder$mean_normalized_rf), ], row.names = FALSE)
cat("\nBest method by independent recorder:\n")
print(best_method, row.names = FALSE)
if (exists("calibrated_recorder")) {
  cat("\nAccuracy using the calibrated 0.04 prime-editing rate:\n")
  print(
    calibrated_recorder[order(calibrated_recorder$mean_normalized_rf), ],
    row.names = FALSE
  )
}
