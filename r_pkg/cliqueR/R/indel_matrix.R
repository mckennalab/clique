#' Build a cell-by-site indel character matrix
#'
#' Converts the per-read `ce` (called-edit) strings that clique writes into the
#' integer character matrix expected by lineage-tracing tree solvers. Each `ce`
#' value is a `_`-joined list of per-target calls (see [read_longread_lineage()]);
#' reads are grouped into cells, the majority call is taken per site, and each
#' distinct edit string is assigned an integer state. `0` is the uncut/`NONE`
#' reference state and `-1` (by default) denotes missing data.
#'
#' The `ce` fields are split on `_` in YAML target order. Multiple events within
#' one target remain a single compound state (their `&`-joined event string is
#' not split). Ties in the per-cell majority vote are resolved by the first
#' value returned by [base::table()].
#'
#' @param lineage_df A data frame from [read_longread_lineage()] /
#'   [read_10x_lineage()] with a `read_name` column and the event column.
#' @param sites Optional character vector naming the cut sites; defaults to
#'   `site1..siteN` inferred from the event strings.
#' @param cell_col Name of the column identifying cells (e.g. `"e0"` for the
#'   extracted cell barcode). If `NULL`, each read is treated as its own cell.
#' @param event_col Name of the called-edit column (default `"ce"`).
#' @param reference_col Name of the reference column, used only to warn when
#'   cells span multiple references (whose site layouts may differ).
#' @param missing_threshold Cells with more than this fraction of missing sites
#'   are dropped.
#' @param missing_value Integer used to encode missing data (default `-1`).
#' @return A list with:
#'   * `matrix`: integer matrix, rows = cells, cols = sites.
#'   * `state_map`: named list mapping each site to its edit-string -> state-int
#'     dictionary (`NONE` is always `0`).
#' @examples
#' \dontrun{
#' lineage <- read_longread_lineage("sample.consensus.bam", tags = c("ce", "e0"))
#' im <- build_indel_matrix(lineage, cell_col = "e0")
#' dim(im$matrix)            # cells x sites
#' head(im$state_map$site1)  # edit string -> integer state
#' }
#' @export
build_indel_matrix <- function(lineage_df,
                               sites = NULL,
                               cell_col = NULL,
                               event_col = "ce",
                               reference_col = "reference",
                               missing_threshold = 0.5,
                               missing_value = -1L) {
  if (!event_col %in% names(lineage_df)) {
    rlang::abort(sprintf("`%s` column not found in lineage_df.", event_col))
  }

  if (!is.null(reference_col) && reference_col %in% names(lineage_df)) {
    refs <- unique(lineage_df[[reference_col]])
    if (length(refs) > 1) {
      rlang::warn(sprintf(
        "lineage_df spans %d references (%s); site columns are assumed comparable across them.",
        length(refs), paste(refs, collapse = ", ")))
    }
  }

  cell <- if (!is.null(cell_col) && cell_col %in% names(lineage_df)) {
    as.character(lineage_df[[cell_col]])
  } else {
    as.character(lineage_df$read_name)
  }

  ev <- as.character(lineage_df[[event_col]])
  split_ev <- strsplit(ifelse(is.na(ev), "", ev), "_", fixed = TRUE)
  n_sites <- max(lengths(split_ev), 0L)
  if (n_sites == 0L) rlang::abort("No events found to build a matrix from.")
  if (is.null(sites)) sites <- paste0("site", seq_len(n_sites))
  if (length(sites) != n_sites) {
    rlang::abort(sprintf("`sites` has length %d but events imply %d sites.",
                         length(sites), n_sites))
  }

  # Per-read site event strings, padding ragged rows with NA (= missing).
  read_site <- t(vapply(split_ev, function(x) {
    length(x) <- n_sites
    x
  }, character(n_sites)))
  read_site <- matrix(read_site, ncol = n_sites, dimnames = list(NULL, sites))
  # An absent / empty call is missing data, not a distinct empty-string allele.
  read_site[!is.na(read_site) & read_site == ""] <- NA_character_

  # One event string per (cell, site) by majority vote across the cell's reads.
  cells <- unique(cell)
  event_mat <- matrix(NA_character_, nrow = length(cells), ncol = n_sites,
                      dimnames = list(cells, sites))
  for (ci in seq_along(cells)) {
    rows <- which(cell == cells[ci])
    for (sj in seq_len(n_sites)) {
      event_mat[ci, sj] <- .mode_ignore_na(read_site[rows, sj])
    }
  }

  # Map edit strings to integer states per site (NONE = 0, others 1..k sorted).
  state_map <- vector("list", n_sites)
  names(state_map) <- sites
  int_mat <- matrix(missing_value, nrow = length(cells), ncol = n_sites,
                    dimnames = list(cells, sites))
  for (sj in seq_len(n_sites)) {
    col <- event_mat[, sj]
    present <- col[!is.na(col)]
    others <- sort(unique(present[present != "NONE"]))
    others_states <- seq_along(others)
    names(others_states) <- others
    mapping <- c(NONE = 0L, others_states)
    state_map[[sj]] <- mapping
    for (ci in seq_along(cells)) {
      v <- col[ci]
      int_mat[ci, sj] <- if (is.na(v)) missing_value else mapping[[v]]
    }
  }

  miss_frac <- rowMeans(int_mat == missing_value)
  keep <- miss_frac <= missing_threshold
  int_mat <- int_mat[keep, , drop = FALSE]

  list(matrix = int_mat, state_map = state_map)
}

#' Most common non-missing value; `NA` if all missing.
#' @noRd
.mode_ignore_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tab <- table(x)
  names(tab)[which.max(tab)]
}

#' Convert an indel character matrix to a normalized Hamming distance matrix
#'
#' Distance between two cells is the fraction of their shared non-missing sites
#' that carry different states (optionally site-weighted). Pairs with no shared
#' non-missing site are assigned distance 1.
#'
#' @param indel_mat Output `$matrix` from [build_indel_matrix()].
#' @param weights Optional numeric vector of per-site weights (length = sites).
#' @param missing Integer value used to encode missing data in `indel_mat`.
#' @return A `dist` object suitable for [ape::nj()].
#' @examples
#' \dontrun{
#' d <- indel_distance(im$matrix)
#' tree <- ape::nj(d)
#' }
#' @export
indel_distance <- function(indel_mat,
                           weights = NULL,
                           missing = -1L) {
  n <- nrow(indel_mat)
  p <- ncol(indel_mat)
  if (is.null(weights)) weights <- rep(1, p)
  if (length(weights) != p) {
    rlang::abort(sprintf("`weights` length %d must equal number of sites %d.",
                         length(weights), p))
  }

  d <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (i == n) break
    for (j in (i + 1):n) {
      vi <- indel_mat[i, ]
      vj <- indel_mat[j, ]
      valid <- vi != missing & vj != missing
      wsum <- sum(weights[valid])
      dij <- if (wsum == 0) 1 else sum(weights[valid] * (vi[valid] != vj[valid])) / wsum
      d[i, j] <- dij
      d[j, i] <- dij
    }
  }
  rownames(d) <- rownames(indel_mat)
  stats::as.dist(d)
}
