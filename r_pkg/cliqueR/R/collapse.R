#' Collapse aligned reads into UMI consensus sequences
#'
#' Wraps `clique collapse`. Takes an aligned BAM and produces a per-UMI
#' consensus BAM via clustering and stretcher-based consensus building.
#'
#' @details
#' This typed wrapper is reserved for the planned high-level API and currently
#' errors without invoking the CLI. Use [generate_clique_script()] with
#' `collapse = TRUE` for a complete workflow, or [clique_run()] to invoke
#' `collapse` directly.
#'
#' @param input_bam Path to an aligned BAM (typically the output of [align()]).
#' @param read_structure Path to the read-structure YAML.
#' @param output_bam Path to write the collapsed/consensus BAM.
#' @param temp_dir Scratch directory for sharded sorting. `NULL` uses the
#'   tool's default.
#' @param max_deletion Maximum deletion length allowed during alignment.
#' @param find_inversions If `TRUE`, also search for inverted matches.
#' @param fast_reference_lookup If `TRUE`, use the kmer-based reference index.
#' @param correct_only If `TRUE`, correct tags without building consensus.
#' @param threads Number of worker threads.
#' @return This development release errors because the wrapper is not yet
#'   implemented.
#' @examples
#' \dontrun{
#' collapse(
#'   input_bam = "sample.aligned.bam",
#'   read_structure = "barcodes.yaml",
#'   output_bam = "sample.consensus.bam",
#'   threads = 8
#' )
#'
#' # Tag correction only — no consensus building
#' collapse(
#'   input_bam = "sample.aligned.bam",
#'   read_structure = "barcodes.yaml",
#'   output_bam = "sample.corrected.bam",
#'   correct_only = TRUE
#' )
#' }
#' @export
collapse <- function(input_bam,
                     read_structure,
                     output_bam,
                     temp_dir = NULL,
                     max_deletion = 0L,
                     find_inversions = FALSE,
                     fast_reference_lookup = FALSE,
                     correct_only = FALSE,
                     threads = 1L) {
  stop("not yet implemented")
}
