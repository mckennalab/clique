#' Read a long-read lineage BAM (PacBio / Nanopore)
#'
#' Parses a BAM produced by the CLI's `clique align` or `clique collapse`
#' command into a per-read data frame, pulling the aux tags clique writes: the
#' called-edit string (`ce`), the read count (`rc`), and any extracted-UMI tags
#' (`e<symbol>`, e.g. `e0` for the cell barcode). Reads are streamed through
#' `samtools view`, so `samtools` must be on the `PATH` (or passed via
#' `samtools =`); this avoids a Bioconductor/Rsamtools dependency.
#'
#' Tag values are returned as character columns regardless of their SAM type.
#' Convert depths such as `rc` and `dc` explicitly before numeric summaries.
#' When `tags = NULL`, tags are discovered from the first BAM record; pass an
#' explicit vector when later records may contain additional annotations.
#'
#' @param bam_path Path to a clique BAM.
#' @param read_structure Optional path to the read-structure YAML used to
#'   produce the BAM. Currently reserved for mapping UMI symbols to tag names;
#'   tag discovery presently reads the BAM directly.
#' @param tags Character vector of BAM tag names to extract. If `NULL`, all aux
#'   tags found on the first record are used (plus `ce` is always included).
#' @param region Optional region string (e.g. `"left"`) passed to `samtools
#'   view` to restrict to one reference; requires an index for that BAM.
#' @param samtools Path to the `samtools` executable.
#' @return A data frame keyed by `read_name` with a `reference` column and one
#'   column per extracted tag (missing tags are `NA`).
#' @seealso `vignette("bam-output", package = "cliqueR")`,
#'   [build_indel_matrix()]
#' @examples
#' \dontrun{
#' reads <- read_longread_lineage("sample.consensus.bam", tags = c("ce", "e0"))
#' }
#' @export
read_longread_lineage <- function(bam_path,
                                  read_structure = NULL,
                                  tags = NULL,
                                  region = NULL,
                                  samtools = "samtools") {
  if (!file.exists(bam_path)) {
    rlang::abort(sprintf("BAM file not found: %s", bam_path))
  }

  argv <- c("view", bam_path)
  if (!is.null(region)) argv <- c(argv, region)

  result <- processx::run(
    command = samtools,
    args = argv,
    error_on_status = FALSE
  )
  if (!identical(result$status, 0L)) {
    rlang::abort(c(
      sprintf("`samtools view` failed (exit %s) on %s", format(result$status), bam_path),
      i = utils::tail(strsplit(result$stderr %||% "", "\n", fixed = TRUE)[[1]], 10)
    ))
  }

  lines <- strsplit(result$stdout %||% "", "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(lines)]

  if (length(lines) == 0) {
    empty <- data.frame(read_name = character(), reference = character(),
                        stringsAsFactors = FALSE)
    for (tg in tags) empty[[tg]] <- character()
    return(empty)
  }

  fields <- strsplit(lines, "\t", fixed = TRUE)
  read_name <- vapply(fields, function(f) f[[1]], character(1))
  reference <- vapply(fields, function(f) f[[3]], character(1))

  # Auto-discover tags from the first record if not specified.
  if (is.null(tags)) {
    opt <- fields[[1]]
    opt <- if (length(opt) > 11) opt[12:length(opt)] else character()
    tags <- sub(":[AifZHBc]:.*$", "", opt)
    tags <- unique(c(tags, "ce"))
  }

  out <- data.frame(read_name = read_name, reference = reference,
                    stringsAsFactors = FALSE)
  for (tg in tags) {
    out[[tg]] <- .extract_sam_tag(lines, tg)
  }
  out
}

#' Extract one aux tag's value from raw SAM lines (vectorized). Returns `NA` for
#' records lacking the tag. Tag values never contain a tab, so a per-line regex
#' anchored on the `TAG:TYPE:` prefix is sufficient.
#' @noRd
.extract_sam_tag <- function(lines, tag) {
  pat <- sprintf("%s:[AifZHBc]:[^\t]*", tag)
  m <- regexpr(pat, lines)
  raw <- regmatches(lines, m)             # only matching lines
  vals <- rep(NA_character_, length(lines))
  vals[m != -1] <- sub(sprintf("^%s:[AifZHBc]:", tag), "", raw)
  vals
}
