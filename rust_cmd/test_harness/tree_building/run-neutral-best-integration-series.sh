#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
CLIQUE_R_DIR=${CLIQUE_R_DIR:-"$(cd -- "$REPO_ROOT/../r_pkg/cliqueR" && pwd)"}
REMOTE_MITO_DIR=${REMOTE_MITO_DIR:-/Users/aaronmck/Desktop/code/remote_mito_clean}
SOURCE_BENCHMARK=${SOURCE_BENCHMARK:-"$REMOTE_MITO_DIR/output/lineage_benchmark_20260819_173705"}
OUTPUT_ROOT=${OUTPUT_ROOT:-"$SCRIPT_DIR/results/neutral_best_integration_series_n250"}
LEVELS=${LEVELS:-5,10,20,30,40,50}
SYSTEMS=${SYSTEMS:-baseline,wt_crispr,prime,palincode,mitochondrial}
SEEDS=${SEEDS:-1,2,3,4,5}
THRESHOLD=${THRESHOLD:-0.20}
WORKERS=${WORKERS:-4}
TREE_WORKERS=${TREE_WORKERS:-8}

mkdir -p "$OUTPUT_ROOT"
printf '%s\t%s\n' \
  levels "$LEVELS" \
  systems "$SYSTEMS" \
  seeds "$SEEDS" \
  normalized_rf_threshold "$THRESHOLD" \
  stopping_statistic mean_across_six_methods \
  tree_shape neutral \
  sampled_cells 250 \
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
  --workers="$WORKERS"

IFS=',' read -r -a level_values <<<"$LEVELS"
active_systems=$SYSTEMS
for level in "${level_values[@]}"; do
  if [[ -z "$active_systems" ]]; then
    echo "All recorders crossed the stopping threshold before level $level."
    break
  fi

  reconstruction_dir="$OUTPUT_ROOT/reconstruction/level_$(printf '%02d' "$level")"
  echo "Reconstructing level $level for: $active_systems"
  Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
    --benchmark-dir="$OUTPUT_ROOT/simulation" \
    --package-dir="$CLIQUE_R_DIR" \
    --output-dir="$reconstruction_dir" \
    --methods=nj,parsimony,iqtree,mix,vine \
    --matrix=character_matrix \
    --shapes=neutral \
    --systems="$active_systems" \
    --integrations="$level" \
    --observation-depths="$level" \
    --seeds="$SEEDS" \
    --workers="$TREE_WORKERS" \
    --threads=1 \
    --vine-nj-only=false \
    --resume=true

  MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/cliqueR-mpl-cache} \
  Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
    --benchmark-dir="$OUTPUT_ROOT/simulation" \
    --package-dir="$CLIQUE_R_DIR" \
    --output-dir="$reconstruction_dir" \
    --methods=nj,parsimony,iqtree,mix,vine,cassiopeia \
    --matrix=character_matrix \
    --shapes=neutral \
    --systems="$active_systems" \
    --integrations="$level" \
    --observation-depths="$level" \
    --seeds="$SEEDS" \
    --workers=1 \
    --threads=1 \
    --vine-nj-only=false \
    --resume=true

  Rscript "$SCRIPT_DIR/summarize-neutral-best-integration-level.R" \
    "$OUTPUT_ROOT" "$level" "$THRESHOLD"
  active_systems=$(tr -d '\r\n' <"$OUTPUT_ROOT/active_systems.txt")
done

Rscript "$SCRIPT_DIR/finalize-neutral-best-integration-series.R" "$OUTPUT_ROOT"
Rscript "$SCRIPT_DIR/plot-neutral-best-integration-series.R" "$OUTPUT_ROOT"
write_complete="$OUTPUT_ROOT/benchmark_complete.txt"
printf 'complete\n' >"$write_complete"
echo "Neutral best-parameter integration series complete: $OUTPUT_ROOT"
