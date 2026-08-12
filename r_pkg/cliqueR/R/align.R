#' Align FASTQ reads to a reference with clique
#'
#' Wraps `clique align`. Takes raw FASTQ inputs plus a read-structure YAML and
#' writes an aligned BAM.
#'
#' @details
#' This typed wrapper is reserved for the planned high-level API and currently
#' errors without invoking the CLI. Use [generate_clique_script()] for a
#' complete align/sort/collapse workflow, or [clique_run()] to invoke `align`
#' directly.
#'
#' @param read1 Path to R1 FASTQ (required).
#' @param read2 Path to R2 FASTQ, or `NULL` for single-end.
#' @param index1,index2 Optional index-read FASTQs.
#' @param read_structure Path to the read-structure YAML.
#' @param output_bam Path to write the aligned BAM.
#' @param aligner One of `"WFA"`, `"Degenerate"`, `"Inversion"`.
#' @param max_reference_multiplier Cap on read-to-reference length ratio.
#' @param min_read_length Reads shorter than this are dropped.
#' @param threads Number of worker threads.
#' @return This development release errors because the wrapper is not yet
#'   implemented.
#' @examples
#' \dontrun{
#' # Paired-end short-read run
#' align(
#'   read1 = "sample_R1.fastq.gz",
#'   read2 = "sample_R2.fastq.gz",
#'   read_structure = "barcodes.yaml",
#'   output_bam = "sample.aligned.bam",
#'   threads = 8
#' )
#'
#' # Long-read PacBio CCS
#' align(
#'   read1 = "sample.ccs.fastq.gz",
#'   read_structure = "longread.yaml",
#'   output_bam = "sample.aligned.bam",
#'   aligner = "WFA"
#' )
#' }
#' @export
align <- function(read1,
                  read2 = NULL,
                  index1 = NULL,
                  index2 = NULL,
                  read_structure,
                  output_bam,
                  aligner = c("WFA", "Degenerate", "Inversion"),
                  max_reference_multiplier = 2L,
                  min_read_length = 50L,
                  threads = 1L) {
  stop("not yet implemented")
}
