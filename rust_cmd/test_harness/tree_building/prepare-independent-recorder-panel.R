#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}

source_root <- normalizePath(value_after(
  "--benchmark-dir",
  "/Users/aaronmck/Desktop/code/remote_mito_clean/output/lineage_benchmark_20260819_173705"
), mustWork = TRUE)
output_root <- normalizePath(value_after(
  "--output-dir",
  "test_harness/tree_building/results/independent_recorders_n250_level05/simulation"
), mustWork = FALSE)
tree_size <- as.integer(value_after("--tree-size", "250"))
level <- as.integer(value_after("--level", "5"))
seed <- as.integer(value_after("--seed", "1"))
shapes <- strsplit(value_after(
  "--shapes", "balanced,comb,neutral,hierarchical"
), ",", fixed = TRUE)[[1L]]
systems <- c("baseline", "prime", "palincode", "mitochondrial")

if (is.na(tree_size) || tree_size < 2L || is.na(level) || level < 1L ||
    is.na(seed) || seed < 0L) {
  stop("Tree size, level, and seed must be valid non-negative integers.")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("The Matrix package is required.")
}

manifest <- read.csv(
  gzfile(file.path(source_root, "benchmark_manifest.csv.gz")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
keep <- manifest$population_seed == seed &
  manifest$tree_size == tree_size &
  manifest$shape %in% shapes &
  manifest$system %in% systems &
  (
    (manifest$system != "mitochondrial" & manifest$integrations == level) |
    (manifest$system == "mitochondrial" & manifest$observation_depth == level)
  )
selected <- manifest[keep, , drop = FALSE]
selected <- selected[order(selected$shape, selected$system), , drop = FALSE]
expected <- length(shapes) * length(systems)
if (nrow(selected) != expected || anyDuplicated(selected[, c("shape", "system")])) {
  stop("Expected exactly one condition for each requested shape and recorder.")
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)
diagnostics <- vector("list", nrow(selected))

for (index in seq_len(nrow(selected))) {
  condition_id <- selected$condition_id[index]
  condition_dir <- file.path(output_root, "conditions", condition_id)
  dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)
  source_matrix <- readRDS(file.path(
    source_root, selected$character_matrix[index]
  ))
  edited_counts <- Matrix::colSums(source_matrix != 0)
  retained <- edited_counts > 0L & edited_counts < nrow(source_matrix)
  filtered_matrix <- source_matrix[, retained, drop = FALSE]
  if (!ncol(filtered_matrix)) {
    stop("Condition has no variable characters: ", condition_id)
  }

  matrix_path <- file.path(condition_dir, "character_matrix_variable.rds")
  truth_path <- file.path(condition_dir, "ground_truth_tree.nwk")
  saveRDS(filtered_matrix, matrix_path)
  if (!file.copy(
    file.path(source_root, selected$ground_truth_tree[index]),
    truth_path,
    overwrite = TRUE
  )) {
    stop("Could not copy the ground-truth tree for ", condition_id)
  }

  selected$character_matrix[index] <- file.path(
    "conditions", condition_id, basename(matrix_path)
  )
  selected$ground_truth_tree[index] <- file.path(
    "conditions", condition_id, basename(truth_path)
  )
  selected$character_matrix_columns[index] <- ncol(filtered_matrix)
  selected$character_matrix_nonzero[index] <- length(filtered_matrix@x)
  diagnostics[[index]] <- data.frame(
    condition_id = condition_id,
    shape = selected$shape[index],
    system = selected$system[index],
    cells = nrow(source_matrix),
    original_characters = ncol(source_matrix),
    retained_variable_characters = ncol(filtered_matrix),
    informative_characters = sum(
      edited_counts[retained] >= 2L &
        edited_counts[retained] <= nrow(source_matrix) - 2L
    ),
    singleton_characters = sum(edited_counts[retained] == 1L),
    unique_profiles = nrow(unique(as.matrix(filtered_matrix))),
    nonzero_fraction = length(filtered_matrix@x) /
      (nrow(filtered_matrix) * ncol(filtered_matrix))
  )
}

manifest_connection <- gzfile(
  file.path(output_root, "benchmark_manifest.csv.gz"), open = "wt"
)
write.csv(selected, manifest_connection, row.names = FALSE)
close(manifest_connection)
diagnostics <- do.call(rbind, diagnostics)
write.table(
  diagnostics,
  file.path(output_root, "matrix_diagnostics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
settings <- list(
  source_benchmark = source_root,
  seed = seed,
  tree_size = tree_size,
  integration_or_recovery_level = level,
  shapes = shapes,
  systems = systems,
  matrix = "character_matrix",
  filtering = "Removed characters constant within each 250-cell condition"
)
saveRDS(settings, file.path(output_root, "benchmark_settings.rds"))
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(
    settings, file.path(output_root, "benchmark_settings.json"),
    pretty = TRUE, auto_unbox = TRUE
  )
}
writeLines("complete", file.path(output_root, "benchmark_complete.txt"))
cat(sprintf("Prepared %d independent recorder conditions under %s\n", nrow(selected), output_root))
