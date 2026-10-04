#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
CLIQUE_R_DIR=${CLIQUE_R_DIR:-"$(cd -- "$REPO_ROOT/../r_pkg/cliqueR" && pwd)"}
REMOTE_MITO_DIR=${REMOTE_MITO_DIR:-/Users/aaronmck/Desktop/code/remote_mito_clean}
SOURCE_BENCHMARK=${SOURCE_BENCHMARK:-"$REMOTE_MITO_DIR/output/lineage_benchmark_20260819_173705"}
OUTPUT_ROOT=${OUTPUT_ROOT:-"$SCRIPT_DIR/results/neutral_best_nj_integration_series_nested_n250"}
LEVELS=${LEVELS:-1,5,10,20,30,40,50,100}
SYSTEMS=${SYSTEMS:-baseline,wt_crispr,prime,palincode,mitochondrial}
SEEDS=${SEEDS:-1,2,3,4,5}
THRESHOLD=${THRESHOLD:-0.20}
MITOCHONDRIAL_GENOMES=${MITOCHONDRIAL_GENOMES:-100}
WORKERS=${WORKERS:-4}
TREE_WORKERS=${TREE_WORKERS:-8}

mkdir -p "$OUTPUT_ROOT"
printf '%s\t%s\n' \
  levels "$LEVELS" \
  systems "$SYSTEMS" \
  seeds "$SEEDS" \
  normalized_rf_threshold "$THRESHOLD" \
  stopping_statistic nj_mean_across_seeds \
  stop_early false \
  integration_nesting maximum_level_prefix \
  mitochondrial_nesting fixed_sampling_order_prefix \
  tree_shape neutral \
  sampled_cells 250 \
  mitochondrial_genomes_per_cell "$MITOCHONDRIAL_GENOMES" \
  >"$OUTPUT_ROOT/run_settings.tsv"

Rscript "$SCRIPT_DIR/generate-recorder-parameter-grid.R" \
  --source-root="$REMOTE_MITO_DIR" \
  --benchmark-root="$SOURCE_BENCHMARK" \
  --output-dir="$OUTPUT_ROOT/simulation" \
  --grid=best \
  --systems="$SYSTEMS" \
  --shapes=neutral \
  --seeds="$SEEDS" \
  --levels="$LEVELS" \
  --mitochondrial-genomes="$MITOCHONDRIAL_GENOMES" \
  --workers="$WORKERS"

IFS=',' read -r -a level_values <<<"$LEVELS"
for level in "${level_values[@]}"; do
  reconstruction_dir="$OUTPUT_ROOT/reconstruction/level_$(printf '%02d' "$level")"
  echo "Reconstructing NJ at level $level for: $SYSTEMS"
  Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
    --benchmark-dir="$OUTPUT_ROOT/simulation" \
    --package-dir="$CLIQUE_R_DIR" \
    --output-dir="$reconstruction_dir" \
    --methods=nj \
    --matrix=character_matrix \
    --shapes=neutral \
    --systems="$SYSTEMS" \
    --integrations="$level" \
    --observation-depths="$level" \
    --seeds="$SEEDS" \
    --workers="$TREE_WORKERS" \
    --threads=1 \
    --resume=true

  Rscript "$SCRIPT_DIR/summarize-neutral-best-integration-level.R" \
    "$OUTPUT_ROOT" "$level" "$THRESHOLD" nj
done

Rscript "$SCRIPT_DIR/finalize-neutral-best-integration-series.R" "$OUTPUT_ROOT"
Rscript "$SCRIPT_DIR/plot-neutral-best-nj-integration-series.R" "$OUTPUT_ROOT"
printf 'complete\n' >"$OUTPUT_ROOT/benchmark_complete.txt"
echo "Neutral best-parameter NJ integration series complete: $OUTPUT_ROOT"
