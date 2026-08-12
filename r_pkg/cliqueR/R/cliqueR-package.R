#' cliqueR: Lineage tracing analysis with the clique toolkit
#'
#' @description
#' `cliqueR` connects the `clique` Rust command-line application to R. The CLI
#' aligns amplicon reads, extracts and corrects capture tags, collapses
#' molecule groups, and annotates CRISPR events. The R package provides
#' configuration and execution helpers plus BAM import, character-matrix,
#' tree-reconstruction, and visualization tools.
#'
#' The Rust binary is installed separately. Use [clique_binary()] to inspect
#' the executable selected for the current session or [set_clique_binary()] to
#' override it.
#'
#' @section Running clique:
#' [fasta_targets_to_yaml()] creates a basic read-structure YAML from reference
#' and target sequences. [generate_clique_script()] writes a reproducible
#' align/sort/collapse shell script. [clique_run()] is the lower-level interface
#' for CLI subcommands and flags not represented by the script generator.
#'
#' The typed [align()] and [collapse()] functions are reserved APIs and are not
#' implemented in this development release.
#'
#' @section Reading and analyzing results:
#' [read_longread_lineage()] imports selected auxiliary tags from an aligned or
#' collapsed BAM. [build_indel_matrix()] turns the `ce` edit tag into integer
#' character states, [indel_distance()] computes pairwise distances while
#' ignoring missing calls, and the `tree_*()` functions reconstruct lineage
#' trees. [tree_nj()] and [tree_parsimony()] run in R; [tree_iqtree()],
#' [tree_vine()], [tree_mix()], and [tree_cassiopeia()] connect optional
#' external backends. Use [plot_lineage_tree()] to add categorical or
#' continuous tip annotations.
#'
#' @section BAM tags:
#' Capture tags are named `e0` through `e9` from their YAML symbols; collapsed
#' BAMs can also contain original values `o0` through `o9`. `ce` stores
#' chemistry-filtered event calls, `pe` stores prime-edit classifications, `rc`
#' is corrected-group depth, `dc` is consensus depth after downsampling, and
#' `ar` lists the selected source-read names. See
#' `vignette("bam-output", package = "cliqueR")` for the complete contract.
#'
#' @section Optional dependencies:
#' Most analysis and visualization dependencies are in `Suggests` and are
#' checked only by the functions that need them. `samtools` must be available
#' to [read_longread_lineage()]. The `clique` Rust executable is a separate
#' system requirement. IQ-TREE 2, VINE, PHYLIP MIX, and the Python Cassiopeia
#' package are optional and only needed by their corresponding tree functions.
#'
#' @seealso
#' [clique_binary()], [generate_clique_script()], [read_longread_lineage()],
#' [build_indel_matrix()], [tree_nj()], [tree_iqtree()], [tree_vine()]
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
