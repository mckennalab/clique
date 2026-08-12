#' Read a 10x CellRanger lineage-tracing output
#'
#' Loads a CellRanger output directory and pairs cell barcodes with their
#' lineage-cassette reads. Intended as input to [build_indel_matrix()].
#'
#' @details
#' This function is a reserved API in the current development release and
#' errors without reading either input. For clique BAM files that already carry
#' extracted cell-barcode tags, use [read_longread_lineage()] and pass the
#' relevant tag (for example, `tags = c("e0", "ce")`).
#'
#' @param cellranger_dir Path to a CellRanger `outs/` directory (must contain
#'   `filtered_feature_bc_matrix/` for cell barcodes).
#' @param lineage_bam Path to the per-cell lineage BAM, typically the output
#'   of [collapse()].
#' @param min_reads_per_cell Drop cells with fewer than this many lineage reads.
#' @return This development release errors because 10x import is not yet
#'   implemented.
#' @examples
#' \dontrun{
#' lineage <- read_10x_lineage(
#'   cellranger_dir = "sample/outs",
#'   lineage_bam    = "sample.consensus.bam",
#'   min_reads_per_cell = 3
#' )
#' }
#' @export
read_10x_lineage <- function(cellranger_dir,
                             lineage_bam,
                             min_reads_per_cell = 1L) {
  stop("not yet implemented")
}
