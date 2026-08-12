# cliqueR 0.0.0.9000 (development)

## Documentation

* Added an end-to-end README and getting-started vignette covering binary
  discovery, read-structure generation, pipeline scripts, BAM import, character
  matrices, neighbor-joining trees, and visualization.
* Added a BAM output vignette documenting capture tags, collapse depth and
  provenance, CRISPR/prime-edit calls, reference-routing annotations, and SAM
  conventions.
* Clearly distinguish implemented functions from reserved APIs that currently
  error.
* Expanded package and function help, and added generated help pages for
  `fasta_targets_to_yaml()`, `generate_clique_script()`, and
  `clique_categorical_palette()`.

## Implemented

* `clique_binary()` / `set_clique_binary()` locate or override the CLI.
* `clique_run()` invokes arbitrary subcommands and flags.
* `fasta_targets_to_yaml()` builds basic read structures from FASTA targets.
* `generate_clique_script()` writes align/sort/collapse pipeline scripts.
* `read_longread_lineage()` imports selected clique BAM tags through `samtools`.
* `build_indel_matrix()` / `indel_distance()` prepare lineage characters.
* `tree_nj()` reconstructs neighbor-joining trees.
* `tree_parsimony()` reconstructs unordered maximum-parsimony trees in R with
  phangorn's parsimony ratchet.
* `tree_iqtree()` runs IQ-TREE 2 maximum-likelihood inference on recoded
  morphological characters.
* `tree_vine()` runs VINE's CRISPR variational-inference model and can retain
  posterior tree samples.
* `tree_cassiopeia()` exposes greedy, ILP, hybrid, UPGMA, and neighbor-joining
  Cassiopeia solvers through reticulate.
* `tree_mix()` runs PHYLIP MIX with one-hot edit characters under Camin-Sokal
  or Wagner parsimony.
* `plot_lineage_tree()` visualizes trees and per-tip annotations.

## Reserved APIs

The typed `align()` / `collapse()` wrappers, 10x import, and tree-comparison
metrics are not yet implemented.
