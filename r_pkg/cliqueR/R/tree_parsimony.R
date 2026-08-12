#' Parsimony tree from an indel character matrix
#'
#' Builds an unordered maximum-parsimony tree from the indel character matrix
#' using phangorn's parsimony ratchet. Unlike [tree_mix()], this implementation
#' runs entirely in R and treats each integer at a site as one categorical
#' state.
#'
#' @param indel_mat Output `$matrix` from [build_indel_matrix()].
#' @param start An optional starting tree (defaults to NJ on Hamming distance).
#' @param trace If `TRUE`, print phangorn search progress.
#' @return An `ape::phylo` object with original cell names as tip labels.
#' @examples
#' \dontrun{
#' tree <- tree_parsimony(im$matrix)
#' }
#' @export
tree_parsimony <- function(indel_mat, start = NULL, trace = FALSE) {
  mat <- .validate_indel_matrix(indel_mat, missing = -1L)
  .require_tree_package("ape", "tree_parsimony")
  .require_tree_package("phangorn", "tree_parsimony")
  if (!is.logical(trace) || length(trace) != 1L || is.na(trace)) {
    rlang::abort("`trace` must be TRUE or FALSE.")
  }

  if (is.null(start)) {
    start <- tree_nj(indel_distance(mat, missing = -1L))
  } else {
    if (!inherits(start, "phylo")) {
      rlang::abort("`start` must be an `ape::phylo` object.")
    }
    if (!setequal(start$tip.label, rownames(mat))) {
      rlang::abort("`start` and `indel_mat` must contain the same cell names.")
    }
  }

  states <- sort(unique(mat[mat != -1L]))
  chars <- matrix(as.character(mat), nrow(mat), ncol(mat), dimnames = dimnames(mat))
  chars[mat == -1L] <- "?"
  data <- phangorn::phyDat(chars, type = "USER", levels = as.character(states))
  tree <- phangorn::pratchet(
    data,
    start = start,
    trace = as.integer(trace)
  )
  .annotate_backend_tree(tree, "phangorn-parsimony")
}
