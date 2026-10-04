#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}

source_root <- normalizePath(value_after(
  "--source-root",
  "/Users/aaronmck/Desktop/code/remote_mito_clean"
), mustWork = TRUE)
benchmark_root <- normalizePath(value_after(
  "--benchmark-root",
  file.path(source_root, "output", "lineage_benchmark_20260819_173705")
), mustWork = TRUE)
output_root <- normalizePath(value_after(
  "--output-dir",
  "test_harness/tree_building/results/prime_rate_grid_n250_k05/simulation"
), mustWork = FALSE)
rates <- as.numeric(strsplit(value_after("--rates", "0.01,0.02,0.04,0.08"), ",", fixed = TRUE)[[1]])
shapes <- strsplit(value_after(
  "--shapes", "balanced,comb,neutral,hierarchical"
), ",", fixed = TRUE)[[1]]
tree_size <- as.integer(value_after("--tree-size", "250"))
integrations <- as.integer(value_after("--integrations", "5"))
population_seed <- as.integer(value_after("--seed", "1"))

if (any(!is.finite(rates)) || any(rates <= 0) || any(rates >= 1)) {
  stop("--rates must contain probabilities strictly between zero and one.")
}
if (is.na(tree_size) || tree_size < 2L || is.na(integrations) || integrations < 1L) {
  stop("--tree-size must be at least 2 and --integrations must be positive.")
}

# The simulator moved into a package layout: its core files now live in
# origin/R and the benchmarking library in analysis/. load_origin.R sources
# both into the global environment, which this script depends on because it
# replaces several simulator functions with its own versions below.
source(file.path(source_root, "load_origin.R"))

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)
base_recording_params <- lineage_benchmark_recording_params
active_rate <- NA_real_

# Keep the benchmark implementation intact while exposing its hard-coded prime
# editing probability as the single experimental variable.
lineage_benchmark_recording_params <- function(population_params,
                                               system,
                                               maximum_integrations) {
  params <- base_recording_params(
    population_params, system, maximum_integrations
  )
  if (identical(tolower(as.character(system)), "prime")) {
    params$prime_editing_backend$induced_edit_probability_per_cell_cycle <-
      active_rate
  }
  params
}

rate_label <- function(rate) sprintf("rate_%03d", as.integer(round(rate * 1000)))
condition_rows <- list()
condition_count <- 0L

for (shape_index in seq_along(shapes)) {
  shape <- shapes[[shape_index]]
  original_run <- file.path(
    benchmark_root,
    paste0("shape_", shape),
    sprintf("seed_%06d", population_seed)
  )
  population_path <- file.path(original_run, "population", "population.rds")
  sample_dir <- file.path(original_run, "samples", paste0("n_", tree_size))
  sample_sets_path <- file.path(original_run, "samples", "sample_sets.rds")
  required <- c(
    population_path,
    sample_sets_path,
    file.path(sample_dir, "ground_truth_tree.nwk"),
    resolve_physicell_csv_path(file.path(sample_dir, "sample_cells.csv"))
  )
  if (any(!file.exists(required))) {
    stop("The source benchmark is missing inputs for shape ", shape, ".")
  }

  population_bundle <- readRDS(population_path)
  sample_sets <- readRDS(sample_sets_path)
  physicell_ids <- sample_sets[[as.character(tree_size)]]
  if (is.null(physicell_ids)) {
    stop("The source benchmark has no n=", tree_size, " sample for ", shape, ".")
  }
  recorder_seed <- as.integer(
    (as.double(population_seed) * 10007 + shape_index * 1009 + 2L * 97) %%
      .Machine$integer.max
  )

  for (rate in rates) {
    active_rate <- rate
    label <- rate_label(rate)
    run_dir <- file.path(
      output_root,
      paste0("shape_", shape),
      sprintf("seed_%06d", population_seed),
      label
    )
    recording_dir <- file.path(
      run_dir, "recorder_prime", sprintf("full_%02d_integrations", integrations)
    )
    if (file.exists(file.path(recording_dir, "recording_complete.txt"))) {
      recording <- load_lineage_benchmark_recorder(recording_dir)
    } else {
      cat(sprintf("Generating %s %s prime recorder at rate %.3f\n", shape, label, rate))
      recording <- simulate_lineage_benchmark_recorder(
        population_bundle = population_bundle,
        system = "prime",
        maximum_integrations = integrations,
        seed = recorder_seed,
        output_dir = recording_dir,
        show_progress = FALSE
      )
    }

    condition_dir <- file.path(
      run_dir,
      "recorder_prime",
      "conditions",
      paste0("n_", tree_size),
      sprintf("k_%02d", integrations)
    )
    condition <- write_lineage_benchmark_condition(
      recording = recording,
      system = "prime",
      shape = shape,
      population_seed = population_seed,
      recorder_seed = recorder_seed,
      physicell_ids = physicell_ids,
      integrations = integrations,
      truth_tree_path = file.path(sample_dir, "ground_truth_tree.nwk"),
      sample_cells_path = file.path(sample_dir, "sample_cells.csv"),
      full_recording_dir = recording_dir,
      output_dir = condition_dir,
      benchmark_root = output_root,
      write_dense_csv = FALSE
    )
    condition$condition_id <- sub(
      "_prime_", paste0("_prime_", label, "_"), condition$condition_id, fixed = TRUE
    )
    condition$editing_rate <- rate
    write_physicell_csv(
      condition,
      file.path(condition_dir, "condition_manifest.csv"),
      row.names = FALSE,
      compress = TRUE
    )
    condition_count <- condition_count + 1L
    condition_rows[[condition_count]] <- condition
  }
}

manifest <- do.call(rbind, condition_rows)
rownames(manifest) <- NULL
write_physicell_csv(
  manifest,
  file.path(output_root, "benchmark_manifest.csv"),
  row.names = FALSE,
  compress = TRUE
)
settings <- list(
  source_benchmark = benchmark_root,
  rates = rates,
  shapes = shapes,
  tree_size = tree_size,
  integrations = integrations,
  population_seed = population_seed,
  paired_recorder_seeds = TRUE
)
saveRDS(settings, file.path(output_root, "benchmark_settings.rds"))
lineage_benchmark_write_json(settings, file.path(output_root, "benchmark_settings.json"))
writeLines("complete", file.path(output_root, "benchmark_complete.txt"))
cat(sprintf("Wrote %d rate-grid conditions under %s\n", nrow(manifest), output_root))
