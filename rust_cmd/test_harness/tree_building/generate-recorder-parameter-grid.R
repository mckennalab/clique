#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
script_argument <- commandArgs(trailingOnly = FALSE)
script_argument <- script_argument[startsWith(script_argument, "--file=")]
if (!length(script_argument)) stop("Unable to resolve this script's directory.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_argument[[1L]])))
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
  "test_harness/tree_building/results/recorder_parameter_grid_n250_level05/simulation"
), mustWork = FALSE)
tree_size <- as.integer(value_after("--tree-size", "250"))
level_values <- strsplit(
  value_after("--levels", value_after("--level", "5")), ",", fixed = TRUE
)[[1L]]
level_values <- trimws(level_values)
levels <- suppressWarnings(as.integer(level_values))
level <- if (length(levels) && !anyNA(levels)) max(levels) else NA_integer_
seed_values <- strsplit(
  value_after("--seeds", value_after("--seed", "1")), ",", fixed = TRUE
)[[1L]]
seed_values <- trimws(seed_values)
population_seeds <- suppressWarnings(as.integer(seed_values))
mitochondrial_genomes <- as.integer(value_after("--mitochondrial-genomes", "50"))
workers <- as.integer(value_after("--workers", "4"))
grid_mode <- tolower(value_after("--grid", "wide"))
# Number of doublings added above each extendable recorder's very_high tier.
extension_doublings <- as.integer(value_after("--extension-doublings", "5"))
# Optional explicit tier selection, applied after the grid mode filter. The
# bracket-search driver uses it to add one tier at a time.
requested_settings <- trimws(strsplit(
  value_after("--settings", ""), ",", fixed = TRUE
)[[1L]])
requested_settings <- requested_settings[nzchar(requested_settings)]
systems <- strsplit(value_after(
  "--systems", "baseline,wt_crispr,prime,palincode,mitochondrial"
), ",", fixed = TRUE)[[1L]]
shapes <- strsplit(value_after(
  "--shapes", "balanced,comb,neutral,hierarchical"
), ",", fixed = TRUE)[[1L]]
canonical_shapes <- c("balanced", "comb", "neutral", "hierarchical")
canonical_systems <- c(
  "baseline", "wt_crispr", "prime", "palincode", "mitochondrial"
)

if (is.na(tree_size) || tree_size < 2L || is.na(level) || level < 1L ||
    !length(levels) || anyNA(levels) || any(levels < 1L) ||
    anyDuplicated(levels) || any(as.character(levels) != level_values) ||
    is.unsorted(levels, strictly = TRUE) ||
    !length(population_seeds) || anyNA(population_seeds) ||
    any(population_seeds < 0L) || anyDuplicated(population_seeds) ||
    any(as.character(population_seeds) != seed_values) ||
    is.na(mitochondrial_genomes) || mitochondrial_genomes < level ||
    is.na(workers) || workers < 1L ||
    any(!systems %in% canonical_systems) || anyDuplicated(systems) ||
    any(!shapes %in% canonical_shapes) || anyDuplicated(shapes) ||
    is.na(extension_doublings) || extension_doublings < 0L ||
    extension_doublings > 12L ||
    !grid_mode %in% c("focused", "wide", "best", "bracket", "extension")) {
  stop(paste(
    "Tree size, level, seed, worker count, mitochondrial genomes,",
    "extension doublings, systems, or shapes are invalid."
  ))
}

# The simulator moved into a package layout: its core files now live in
# origin/R and the benchmarking library in analysis/. load_origin.R sources
# both into the global environment, which this script depends on because it
# replaces several simulator functions with its own versions below.
source(file.path(source_root, "load_origin.R"))
source(file.path(script_dir, "flare-wt-crispr-recorder.R"))
# Sourcing the recorder registers it. The simulator consults the registry at
# each hook point, so the four wrappers that used to be assigned over
# prepare_physicell_recording_model(), mutate_physicell_barcode_segment(),
# physicell_baseline_target_layout() and lineage_benchmark_character_matrix()
# here are gone: there is nothing left to forget to install.
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("The Matrix package is required to generate recorder conditions.")
}

parameter_grid <- list(
  list(
    system = "baseline", setting = "very_low",
    parameter = "be_mutations_per_target_per_division",
    primary = 0.00238, secondary = NA_real_
  ),
  list(
    system = "baseline", setting = "low",
    parameter = "be_mutations_per_target_per_division",
    primary = 0.00476, secondary = NA_real_
  ),
  list(
    system = "baseline", setting = "mid",
    parameter = "be_mutations_per_target_per_division",
    primary = 0.00952, secondary = NA_real_
  ),
  list(
    system = "baseline", setting = "high",
    parameter = "be_mutations_per_target_per_division",
    primary = 0.01904, secondary = NA_real_
  ),
  list(
    system = "baseline", setting = "very_high",
    parameter = "be_mutations_per_target_per_division",
    primary = 0.03808, secondary = NA_real_
  ),
  list(
    system = "wt_crispr", setting = "very_low",
    parameter = "nuclease_indel_probabilities_per_target_per_division",
    primary = 0.0025, secondary = 0.025
  ),
  list(
    system = "wt_crispr", setting = "low",
    parameter = "nuclease_indel_probabilities_per_target_per_division",
    primary = 0.005, secondary = 0.05
  ),
  list(
    system = "wt_crispr", setting = "mid",
    parameter = "nuclease_indel_probabilities_per_target_per_division",
    primary = 0.01, secondary = 0.10
  ),
  list(
    system = "wt_crispr", setting = "high",
    parameter = "nuclease_indel_probabilities_per_target_per_division",
    primary = 0.02, secondary = 0.20
  ),
  list(
    system = "wt_crispr", setting = "very_high",
    parameter = "nuclease_indel_probabilities_per_target_per_division",
    primary = 0.04, secondary = 0.40
  ),
  # PEtracer (doi:10.1126/science.adx3800) replaces the single-outcome prime
  # editor in this slot: 3 edit sites per cassette, each installing one of 8
  # predefined marks, at the panel's 5 integrations.
  #
  # The ladder is anchored on edit saturation rather than on the published
  # per-day rate, because the paper's own accuracy criterion is 60-80%
  # saturation, and because the harness records over 7-13 cell cycles against
  # PEtracer's 2-6 week experiment -- transplanting 0.05-0.1 per day directly
  # would sit at 31-53% saturation and understate the recorder.
  #
  # Saturation is 1 - (1 - p)^(efficiency * depth), and the trees are
  # ultrametric at depths 7.00 (neutral), 9.29 (balanced), 9.98 (comb) and
  # 12.88 (hierarchical). The rate reaching 70% varies only 1.74x across
  # shapes, less than one doubling, so one ladder serves all four.
  #
  # The ladder sits one doubling BELOW the 70% anchor, putting the 60-80% band
  # at M-H rather than at M. Centring it exactly on the band was tried and
  # wasted the top tier: 4x the anchor reached 99.6-100% saturation on every
  # shape, and on hierarchical seed 5 it collapsed 250 cells to 2 distinct
  # profiles, which is past the point any method can reconstruct. As set here
  # the band stays interior and every tier stays measurable.
  list(
    system = "prime", setting = "very_low",
    parameter = "induced_edit_probability_per_cell_cycle",
    primary = 0.0225, secondary = NA_real_
  ),
  list(
    system = "prime", setting = "low",
    parameter = "induced_edit_probability_per_cell_cycle",
    primary = 0.045, secondary = NA_real_
  ),
  list(
    system = "prime", setting = "mid",
    parameter = "induced_edit_probability_per_cell_cycle",
    primary = 0.09, secondary = NA_real_
  ),
  list(
    system = "prime", setting = "high",
    parameter = "induced_edit_probability_per_cell_cycle",
    primary = 0.18, secondary = NA_real_
  ),
  list(
    system = "prime", setting = "very_high",
    parameter = "induced_edit_probability_per_cell_cycle",
    primary = 0.36, secondary = NA_real_
  ),
  list(
    system = "palincode", setting = "very_low",
    parameter = "induced_cbit_probabilities_per_cell_cycle",
    primary = 0.02, secondary = 0.0125
  ),
  list(
    system = "palincode", setting = "low",
    parameter = "induced_cbit_probabilities_per_cell_cycle",
    primary = 0.04, secondary = 0.025
  ),
  list(
    system = "palincode", setting = "mid",
    parameter = "induced_cbit_probabilities_per_cell_cycle",
    primary = 0.08, secondary = 0.05
  ),
  list(
    system = "palincode", setting = "high",
    parameter = "induced_cbit_probabilities_per_cell_cycle",
    primary = 0.16, secondary = 0.10
  ),
  list(
    system = "palincode", setting = "very_high",
    parameter = "induced_cbit_probabilities_per_cell_cycle",
    primary = 0.32, secondary = 0.20
  ),
  list(
    system = "mitochondrial", setting = "very_low",
    parameter = "substitution_probability_per_site_per_division",
    primary = 5e-7, secondary = NA_real_
  ),
  list(
    system = "mitochondrial", setting = "low",
    parameter = "substitution_probability_per_site_per_division",
    primary = 1e-6, secondary = NA_real_
  ),
  list(
    system = "mitochondrial", setting = "mid",
    parameter = "substitution_probability_per_site_per_division",
    primary = 2e-6, secondary = NA_real_
  ),
  list(
    system = "mitochondrial", setting = "high",
    parameter = "substitution_probability_per_site_per_division",
    primary = 4e-6, secondary = NA_real_
  ),
  list(
    system = "mitochondrial", setting = "very_high",
    parameter = "substitution_probability_per_site_per_division",
    primary = 8e-6, secondary = NA_real_
  )
)
# Editing-rate extension above very_high.
#
# In the original five-tier grid, BASELINE and mitochondrial recording both
# selected very_high, the top of the tested range, so neither optimum was
# bracketed: the accuracy minimum sat on the edge of the grid and the true best
# rate could lie anywhere above it. These tiers continue each recorder's
# doubling ladder upward so the minimum becomes interior. A tier is named for
# its multiple of that recorder's very_high value, so very_high_4x is four
# times the very_high rate.
#
# The ladder stops before any editing value reaches maximum_editing_probability.
# These values are per-division probabilities, and physicell_probability_hazard()
# silently clamps its input to just under 1, so a tier at or above 1 would
# saturate every target rather than fail. Stopping short keeps every tested
# tier a rate the simulator actually honors.
extension_systems <- c("baseline", "mitochondrial")
maximum_editing_probability <- 0.75
extension_multipliers <- 2^seq_len(extension_doublings)

extension_grid <- list()
for (extension_system in extension_systems) {
  anchor <- Filter(
    function(config) {
      identical(config$system, extension_system) &&
        identical(config$setting, "very_high")
    },
    parameter_grid
  )
  if (!length(anchor)) {
    stop(sprintf(
      "No very_high anchor tier is defined for %s.", extension_system
    ))
  }
  anchor <- anchor[[1L]]
  for (multiplier in extension_multipliers) {
    candidate <- list(
      system = anchor$system,
      setting = sprintf("very_high_%dx", multiplier),
      parameter = anchor$parameter,
      primary = anchor$primary * multiplier,
      secondary = if (is.na(anchor$secondary)) {
        NA_real_
      } else {
        anchor$secondary * multiplier
      }
    )
    candidate_values <- c(candidate$primary, candidate$secondary)
    candidate_values <- candidate_values[!is.na(candidate_values)]
    if (any(candidate_values >= maximum_editing_probability)) break
    extension_grid[[length(extension_grid) + 1L]] <- candidate
  }
}
parameter_grid <- c(parameter_grid, extension_grid)
extension_settings <- unique(vapply(
  extension_grid, function(config) config$setting, character(1)
))
setting_levels <- c(
  "very_low", "low", "mid", "high", "very_high",
  sprintf("very_high_%dx", extension_multipliers)
)

if (grid_mode == "focused") {
  parameter_grid <- Filter(
    function(config) config$setting %in% c("low", "mid", "high"),
    parameter_grid
  )
}
if (grid_mode == "best") {
  best_settings <- c(
    baseline = "mid",
    wt_crispr = "low",
    # PEtracer's optimum moved to high (0.18) in the recorder parameter grid,
    # on the neutral shape this series uses as well as overall.
    prime = "high",
    palincode = "mid",
    mitochondrial = "very_high"
  )
  parameter_grid <- Filter(
    function(config) identical(config$setting, unname(best_settings[[config$system]])),
    parameter_grid
  )
}
# bracket: the complete ladder for the two unbracketed recorders, original
# tiers included, for a self-contained curve.
if (grid_mode == "bracket") {
  parameter_grid <- Filter(
    function(config) config$system %in% extension_systems,
    parameter_grid
  )
}
# extension: only the tiers above very_high, for adding to an existing
# five-tier run without resimulating what it already holds.
if (grid_mode == "extension") {
  parameter_grid <- Filter(
    function(config) {
      config$system %in% extension_systems &&
        config$setting %in% extension_settings
    },
    parameter_grid
  )
}
if (length(requested_settings)) {
  unknown_settings <- setdiff(requested_settings, setting_levels)
  if (length(unknown_settings)) {
    stop(sprintf(
      "Unknown --settings values: %s. Known tiers are: %s.",
      paste(unknown_settings, collapse = ", "),
      paste(setting_levels, collapse = ", ")
    ))
  }
  parameter_grid <- Filter(
    function(config) config$setting %in% requested_settings,
    parameter_grid
  )
}
if (!length(parameter_grid)) {
  stop("The requested grid mode, systems, and settings select no conditions.")
}
parameter_grid <- Filter(
  function(config) config$system %in% systems,
  parameter_grid
)
system_indices <- c(
  baseline = 1L, prime = 2L, palincode = 3L, mitochondrial = 4L,
  wt_crispr = 5L
)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)
base_recording_params <- lineage_benchmark_recording_params
base_mito_params <- lineage_benchmark_mito_params
base_targets_per_integration <- lineage_benchmark_targets_per_integration
base_logical_target_matrix <- lineage_benchmark_logical_target_matrix
base_logical_target_manifest <- lineage_benchmark_logical_target_manifest
base_simulate_recorder <- simulate_lineage_benchmark_recorder
active_config <- NULL

lineage_benchmark_targets_per_integration <- function(system) {
  if (tolower(as.character(system)) == "wt_crispr") {
    10L
  } else {
    base_targets_per_integration(system)
  }
}

lineage_benchmark_logical_target_matrix <- function(state_matrix,
                                                    model,
                                                    integrations) {
  if (isTRUE(model$is_wt_crispr)) {
    flare_wt_crispr_logical_target_matrix(
      state_matrix, model, integrations
    )
  } else {
    base_logical_target_matrix(state_matrix, model, integrations)
  }
}

lineage_benchmark_logical_target_manifest <- function(target_layout, model) {
  if (isTRUE(model$is_wt_crispr)) {
    annotate_flare_wt_crispr_manifest(target_layout)
  } else {
    base_logical_target_manifest(target_layout, model)
  }
}

lineage_benchmark_recording_params <- function(population_params,
                                               system,
                                               maximum_integrations) {
  system <- tolower(as.character(system))
  params <- if (system == "wt_crispr") {
    configure_flare_wt_crispr_params(
      base_recording_params(
        population_params, "baseline", maximum_integrations
      ),
      insertion_probability = active_config$primary,
      deletion_probability = active_config$secondary,
      dropout_probability = 0.30,
      dropout_radius = 200L,
      # Resection fitted to the GSM8791703 HL60 recording: 5.5M deletion events
      # split at the 26bp target spacing into a local class repaired in place
      # (33.5%, mean 8.7bp) and a resecting class (mean 130.5bp, sd 45). Without
      # this the recorder falls back to the legacy mean-one exponential, where a
      # deletion removes a single base and cannot reach the next target -- which
      # is the behaviour the measured lengths rule out.
      deletion_extent = list(
        left = list(shape = 4.1475, rate = 0.06407),
        right = list(shape = 4.1475, rate = 0.06407),
        local = list(probability = 0.335, shape = 0.5931, rate = 0.15473)
      )
    )
  } else {
    base_recording_params(
      population_params, system, maximum_integrations
    )
  }
  if (system == "baseline") {
    for (cell_type in names(params$cell_type_dict$cell_type_params)) {
      params$cell_type_dict$cell_type_params[[cell_type]]$
        induced_editing_params$be_mutations_per_target_per_division <-
        active_config$primary
    }
  } else if (system == "prime") {
    params$prime_editing_backend$induced_edit_probability_per_cell_cycle <-
      active_config$primary
  } else if (system == "palincode") {
    params$palincode_adapter$
      induced_edit_probability_per_cbit_per_cell_cycle <- list(
        active_config$primary, active_config$secondary
      )
  }
  params
}

simulate_lineage_benchmark_recorder <- function(
    population_bundle,
    system,
    maximum_integrations,
    seed,
    output_dir,
    show_progress = TRUE,
    observation_depths = NULL,
    mitochondrial_genomes_per_cell = 32L) {
  recording <- base_simulate_recorder(
    population_bundle = population_bundle,
    system = system,
    maximum_integrations = maximum_integrations,
    seed = seed,
    output_dir = output_dir,
    show_progress = show_progress,
    observation_depths = observation_depths,
    mitochondrial_genomes_per_cell = mitochondrial_genomes_per_cell
  )
  if (tolower(as.character(system)) != "wt_crispr") return(recording)

  manifest <- read_physicell_csv(file.path(output_dir, "run_manifest.csv"))
  manifest$value[manifest$property == "num_targets_per_integration"] <- 10L
  manifest$value[manifest$property == "total_integrated_target_sites"] <-
    10L * maximum_integrations
  write_physicell_csv(
    manifest,
    file.path(output_dir, "run_manifest.csv"),
    row.names = FALSE,
    compress = TRUE
  )
  recording
}

lineage_benchmark_mito_params <- function(population_params,
                                          substitution_probability = 2e-6) {
  base_mito_params(
    population_params,
    substitution_probability = active_config$primary
  )
}

process_seed_shape <- function(task_index) {
  condition_rows <- list()
  diagnostic_rows <- list()
  condition_count <- 0L
  population_seed <- task_grid$population_seed[[task_index]]
  shape <- task_grid$shape[[task_index]]
  shape_seed_index <- match(shape, canonical_shapes)
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

  for (config in parameter_grid) {
    active_config <<- config
    system <- config$system
    recorder_seed <- as.integer(
      (as.double(population_seed) * 10007 + shape_seed_index * 1009 +
         system_indices[[system]] * 97) %% .Machine$integer.max
    )
    run_dir <- file.path(
      output_root,
      paste0("shape_", shape),
      sprintf("seed_%06d", population_seed),
      paste0("recorder_", system),
      paste0("parameter_", config$setting)
    )
    recording_dir <- file.path(run_dir, "full_recording")
    if (file.exists(file.path(recording_dir, "recording_complete.txt"))) {
      recording <- load_lineage_benchmark_recorder(recording_dir)
    } else {
      cat(sprintf(
        "Generating %-12s %-15s %-4s at %g%s\n",
        shape, system, config$setting, config$primary,
        if (is.na(config$secondary)) "" else paste0(",", config$secondary)
      ))
      recording <- simulate_lineage_benchmark_recorder(
        population_bundle = population_bundle,
        system = system,
        maximum_integrations = level,
        seed = recorder_seed,
        output_dir = recording_dir,
        show_progress = FALSE,
        observation_depths = if (system == "mitochondrial") levels else NULL,
        mitochondrial_genomes_per_cell = mitochondrial_genomes
      )
    }

    for (condition_level in levels) {
      condition_dir <- file.path(
        run_dir,
        if (length(levels) == 1L) {
          "condition"
        } else {
          sprintf("condition_level_%02d", condition_level)
        }
      )
      condition <- write_lineage_benchmark_condition(
        recording = recording,
        system = system,
        shape = shape,
        population_seed = population_seed,
        recorder_seed = recorder_seed,
        physicell_ids = physicell_ids,
        integrations = condition_level,
        truth_tree_path = file.path(sample_dir, "ground_truth_tree.nwk"),
        sample_cells_path = file.path(sample_dir, "sample_cells.csv"),
        full_recording_dir = recording_dir,
        output_dir = condition_dir,
        benchmark_root = output_root,
        write_dense_csv = FALSE
      )

      matrix_path <- file.path(
        condition_dir, "recording_character_matrix_sparse.rds"
      )
      matrix_value <- readRDS(matrix_path)
      # A character is informative when the sampled cells do not all agree on
      # it. Counting non-zero entries is only equivalent to that for a BINARY
      # recorder: with a mark alphabet, a site edited in every cell still
      # separates cells whenever they carry different marks, and the old test
      # discarded exactly those columns. For 0/1 matrices this is identical to
      # `0 < edited < n`, so no other recorder changes.
      dense_value <- as.matrix(matrix_value)
      distinct_states <- apply(dense_value, 2L, function(column) {
        length(unique(column))
      })
      retained <- distinct_states > 1L
      # Still reported below as a diagnostic; for a mark alphabet it counts
      # edited cells rather than minor-state cells, so read it as occupancy.
      edited_counts <- Matrix::colSums(matrix_value != 0)
      matrix_value <- matrix_value[, retained, drop = FALSE]
      if (!ncol(matrix_value)) {
        stop(
          "No variable characters for ", shape, " ", system, " ",
          config$setting, " at level ", condition_level
        )
      }
      saveRDS(matrix_value, matrix_path)

      condition$condition_id <- sub(
        paste0("_", system, "_"),
        paste0("_", system, "_param_", config$setting, "_"),
        condition$condition_id,
        fixed = TRUE
      )
      condition$parameter_setting <- config$setting
      condition$parameter_rank <- match(config$setting, setting_levels)
      condition$editing_parameter <- config$parameter
      condition$parameter_value <- config$primary
      condition$parameter_value_secondary <- config$secondary
      condition$character_matrix_columns <- ncol(matrix_value)
      condition$character_matrix_nonzero <- length(matrix_value@x)
      write_physicell_csv(
        condition,
        file.path(condition_dir, "condition_manifest.csv"),
        row.names = FALSE,
        compress = TRUE
      )

      condition_count <- condition_count + 1L
      condition_rows[[condition_count]] <- condition
      diagnostic_rows[[condition_count]] <- data.frame(
        condition_id = condition$condition_id,
        shape = shape,
        system = system,
        parameter_setting = config$setting,
        parameter_rank = match(config$setting, setting_levels),
        parameter_value = config$primary,
        parameter_value_secondary = config$secondary,
        level = condition_level,
        cells = nrow(matrix_value),
        variable_characters = ncol(matrix_value),
        informative_characters = sum(
          edited_counts[retained] >= 2L &
            edited_counts[retained] <= nrow(matrix_value) - 2L
        ),
        singleton_characters = sum(edited_counts[retained] == 1L),
        unique_profiles = nrow(unique(as.matrix(matrix_value))),
        nonzero_fraction = length(matrix_value@x) /
          (nrow(matrix_value) * ncol(matrix_value))
      )
    }
  }
  list(condition_rows = condition_rows, diagnostic_rows = diagnostic_rows)
}

task_grid <- expand.grid(
  population_seed = population_seeds,
  shape = shapes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
task_indices <- seq_len(nrow(task_grid))
if (workers > 1L && length(task_indices) > 1L && .Platform$OS.type != "windows") {
  task_results <- parallel::mclapply(
    task_indices,
    process_seed_shape,
    mc.cores = min(workers, length(task_indices)),
    mc.preschedule = TRUE
  )
} else {
  task_results <- lapply(task_indices, process_seed_shape)
}
failed_tasks <- vapply(task_results, inherits, logical(1), what = "try-error")
if (any(failed_tasks)) {
  labels <- sprintf(
    "%s seed %d",
    task_grid$shape[failed_tasks],
    task_grid$population_seed[failed_tasks]
  )
  stop("Recorder generation failed for: ", paste(labels, collapse = ", "))
}
condition_rows <- unlist(
  lapply(task_results, `[[`, "condition_rows"), recursive = FALSE
)
diagnostic_rows <- unlist(
  lapply(task_results, `[[`, "diagnostic_rows"), recursive = FALSE
)

manifest <- do.call(rbind, condition_rows)
rownames(manifest) <- NULL
manifest_connection <- gzfile(
  file.path(output_root, "benchmark_manifest.csv.gz"), open = "wt"
)
write.csv(manifest, manifest_connection, row.names = FALSE)
close(manifest_connection)
diagnostics <- do.call(rbind, diagnostic_rows)
write.table(
  diagnostics,
  file.path(output_root, "matrix_diagnostics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
settings <- list(
  source_benchmark = benchmark_root,
  population_seeds = population_seeds,
  replicates = length(population_seeds),
  tree_size = tree_size,
  integration_or_recovery_level = if (length(levels) == 1L) level else NULL,
  integration_or_recovery_levels = levels,
  maximum_integration_or_recovery_level = level,
  mitochondrial_genomes_per_cell = mitochondrial_genomes,
  workers = workers,
  grid_mode = grid_mode,
  extension_systems = extension_systems,
  extension_doublings = extension_doublings,
  extension_settings = extension_settings,
  maximum_editing_probability = maximum_editing_probability,
  requested_settings = if (length(requested_settings)) requested_settings else NULL,
  setting_levels = setting_levels,
  systems = systems,
  shapes = shapes,
  parameter_grid = parameter_grid,
  paired_recorder_seeds = TRUE,
  filtering = sprintf(
    paste(
      "Retained characters with nonzero states in at least one but not all",
      "%d sampled cells"
    ),
    tree_size
  )
)
saveRDS(settings, file.path(output_root, "benchmark_settings.rds"))
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(
    settings,
    file.path(output_root, "benchmark_settings.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null"
  )
}
writeLines("complete", file.path(output_root, "benchmark_complete.txt"))
cat(sprintf("Wrote %d recorder-parameter conditions under %s\n", nrow(manifest), output_root))
