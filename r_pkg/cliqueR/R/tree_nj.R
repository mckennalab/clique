#' Neighbor-joining tree from a lineage distance matrix
#'
#' Thin wrapper around [ape::nj()] / [ape::bionj()] that accepts the output of
#' [indel_distance()] and labels tips with cell barcodes.
#'
#' @param d A `dist` from [indel_distance()].
#' @param method One of `"nj"` or `"bionj"`.
#' @return An `ape::phylo` object.
#' @examples
#' \dontrun{
#' d <- indel_distance(im$matrix)
#' tree <- tree_nj(d, method = "bionj")
#' }
#' @export
tree_nj <- function(d, method = c("nj", "bionj")) {
  method <- match.arg(method)
  if (!inherits(d, "dist")) {
    rlang::abort("`d` must be a 'dist' object (e.g. from indel_distance()).")
  }
  if (attr(d, "Size") < 3L) {
    rlang::abort("Neighbor-joining needs at least 3 cells.")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    rlang::abort("Package 'ape' is required for tree_nj(); install it with install.packages('ape').")
  }
  switch(method,
    nj = ape::nj(d),
    bionj = ape::bionj(d)
  )
}
