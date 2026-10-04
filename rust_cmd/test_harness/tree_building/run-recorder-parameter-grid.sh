#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
CLIQUE_R_DIR=${CLIQUE_R_DIR:-"$(cd -- "$REPO_ROOT/../r_pkg/cliqueR" && pwd)"}
REMOTE_MITO_DIR=${REMOTE_MITO_DIR:-/Users/aaronmck/Desktop/code/remote_mito_clean}
SOURCE_BENCHMARK=${SOURCE_BENCHMARK:-"$REMOTE_MITO_DIR/output/lineage_benchmark_20260819_173705"}
OUTPUT_ROOT=${OUTPUT_ROOT:-"$SCRIPT_DIR/results/recorder_parameter_grid_n250_level05"}
WORKERS=${WORKERS:-4}
TREE_WORKERS=${TREE_WORKERS:-8}
GRID=${GRID:-wide}
SYSTEMS=${SYSTEMS:-baseline,wt_crispr,prime,palincode,mitochondrial}
SEEDS=${SEEDS:-1,2,3,4,5}

Rscript "$SCRIPT_DIR/generate-recorder-parameter-grid.R" \
  --source-root="$REMOTE_MITO_DIR" \
  --benchmark-root="$SOURCE_BENCHMARK" \
  --output-dir="$OUTPUT_ROOT/simulation" \
  --grid="$GRID" \
  --systems="$SYSTEMS" \
  --seeds="$SEEDS" \
  --workers="$WORKERS"

Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
  --benchmark-dir="$OUTPUT_ROOT/simulation" \
  --package-dir="$CLIQUE_R_DIR" \
  --output-dir="$OUTPUT_ROOT/reconstruction" \
  --methods=nj,parsimony,iqtree,mix,vine \
  --matrix=character_matrix \
  --seeds="$SEEDS" \
  --workers="$TREE_WORKERS" \
  --threads=1 \
  --vine-nj-only=false \
  --resume=true

# Cassiopeia's Python backend is not stable in forked macOS R workers.
MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/cliqueR-mpl-cache} \
Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
  --benchmark-dir="$OUTPUT_ROOT/simulation" \
  --package-dir="$CLIQUE_R_DIR" \
  --output-dir="$OUTPUT_ROOT/reconstruction" \
  --methods=nj,parsimony,iqtree,mix,vine,cassiopeia \
  --matrix=character_matrix \
  --seeds="$SEEDS" \
  --workers=1 \
  --threads=1 \
  --vine-nj-only=false \
  --resume=true

Rscript "$SCRIPT_DIR/summarize-recorder-parameter-grid.R" "$OUTPUT_ROOT"
Rscript "$SCRIPT_DIR/plot-recorder-parameter-grid.R" "$OUTPUT_ROOT"
