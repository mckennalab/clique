#' Robinson-Foulds distance between two trees
#'
#' This is a reserved API in the current development release and currently
#' errors.
#'
#' @param tree1,tree2 `ape::phylo` objects sharing a tip-label set.
#' @param normalize If `TRUE`, divide by the maximum possible RF (2*(n-3)).
#' @return This development release errors because the metric is not yet
#'   implemented.
#' @examples
#' \dontrun{
#' rf_distance(reconstructed, ground_truth)
#' }
#' @export
rf_distance <- function(tree1, tree2, normalize = TRUE) {
  stop("not yet implemented")
}

#' Triplet correctness between a reconstructed tree and a ground-truth tree
#'
#' Fraction of cell triplets whose resolved topology matches the ground truth.
#'
#' This is a reserved API in the current development release and currently
#' errors.
#'
#' @param reconstructed,ground_truth `ape::phylo` objects.
#' @param n_triplets Number of random triplets to sample. `Inf` enumerates all.
#' @return This development release errors because the metric is not yet
#'   implemented.
#' @examples
#' \dontrun{
#' triplet_correctness(reconstructed, ground_truth, n_triplets = 1e5)
#' }
#' @export
triplet_correctness <- function(reconstructed,
                                ground_truth,
                                n_triplets = 1e4) {
  stop("not yet implemented")
}

#' Quartet distance between two trees
#'
#' This is a reserved API in the current development release and currently
#' errors.
#'
#' @param tree1,tree2 `ape::phylo` objects.
#' @return This development release errors because the metric is not yet
#'   implemented.
#' @export
quartet_distance <- function(tree1, tree2) {
  stop("not yet implemented")
}

#' Ancestor recall: fraction of ground-truth ancestors recovered
#'
#' For each internal node of `ground_truth`, checks whether the same tip-set
#' partition exists in `reconstructed`.
#'
#' This is a reserved API in the current development release and currently
#' errors.
#'
#' @param reconstructed,ground_truth `ape::phylo` objects.
#' @return This development release errors because the metric is not yet
#'   implemented.
#' @export
ancestor_recall <- function(reconstructed, ground_truth) {
  stop("not yet implemented")
}
