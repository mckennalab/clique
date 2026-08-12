#' Build lineage trees with several methods at once
#'
#' Runs any combination of the `cliqueR` tree reconstructions on one indel
#' character matrix and returns the resulting trees, optionally writing a Newick
#' file (`tree_<method>.nwk`) for each. It is a thin dispatcher over
#' [tree_nj()], [tree_parsimony()], [tree_iqtree()], [tree_mix()], [tree_vine()]
#' and [tree_cassiopeia()] — the neighbor-joining distance is built internally
#' from [indel_distance()], every other method receives `indel_mat` directly.
#'
#' Methods whose backend is unavailable (a missing external binary, or
#' Cassiopeia's Python environment) fail per-method: with `on_error = "warn"`
#' (the default) that method is skipped with a warning and returned as `NULL`
#' while the others still run, so a bare call with the two dependency-free
#' methods works everywhere.
#'
#' @param indel_mat A cell x site indel matrix, e.g. from [build_indel_matrix()].
#' @param methods Character vector of methods to run; any of `"nj"`,
#'   `"parsimony"`, `"iqtree"`, `"mix"`, `"vine"`, `"cassiopeia"`. Defaults to the
#'   two methods that need no external dependency.
#' @param output_dir Directory to write one `tree_<method>.nwk` per successful
#'   method (created if needed). `NULL` (default) writes nothing and returns the
#'   trees only.
#' @param args Named list of extra arguments per method, e.g.
#'   `list(nj = list(method = "bionj"), cassiopeia = list(solver = "hybrid"))`.
#' @param missing Missing-state value used for the neighbor-joining distance and
#'   up-front matrix validation. Default `-1L`.
#' @param on_error `"warn"` (skip a failed method, default) or `"stop"` (abort on
#'   the first failure).
#' @return A named list (one entry per requested method) of `ape::phylo` objects
#'   (or `ape::multiPhylo` for methods returning several trees); `NULL` for any
#'   method that was skipped or failed.
#' @examples
#' \dontrun{
#' im <- build_indel_matrix(lineage_df)
#' trees <- build_trees(im$matrix,
#'                      methods = c("nj", "parsimony", "cassiopeia"),
#'                      args = list(cassiopeia = list(solver = "greedy")),
#'                      output_dir = "trees/")
#' rf_distance(trees$nj, trees$parsimony)
#' }
#' @seealso [tree_nj()], [tree_parsimony()], [tree_iqtree()], [tree_mix()],
#'   [tree_vine()], [tree_cassiopeia()], [rf_distance()]
#' @export
build_trees <- function(indel_mat,
                        methods = c("nj", "parsimony"),
                        output_dir = NULL,
                        args = list(),
                        missing = -1L,
                        on_error = c("warn", "stop")) {
  on_error <- match.arg(on_error)
  known <- c("nj", "parsimony", "iqtree", "mix", "vine", "cassiopeia")

  if (!is.character(methods) || length(methods) == 0L) {
    rlang::abort("`methods` must be a non-empty character vector.")
  }
  methods <- unique(methods)
  unknown <- setdiff(methods, known)
  if (length(unknown)) {
    rlang::abort(sprintf(
      "Unknown tree method(s): %s. Known methods: %s.",
      paste(unknown, collapse = ", "), paste(known, collapse = ", ")))
  }
  if (!is.list(args) || (length(args) > 0L && is.null(names(args)))) {
    rlang::abort("`args` must be a named list keyed by method name.")
  }
  bad_args <- setdiff(names(args), known)
  if (length(bad_args)) {
    rlang::abort(sprintf(
      "`args` names must be method names; unexpected: %s.",
      paste(bad_args, collapse = ", ")))
  }
  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1L) {
      rlang::abort("`output_dir` must be a single directory path or NULL.")
    }
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # fail fast on a malformed matrix; reuse the cleaned matrix for every method
  mat <- .validate_indel_matrix(indel_mat, missing = missing)

  build_one <- function(m) {
    extra <- args[[m]]
    if (is.null(extra)) extra <- list()
    switch(m,
      nj         = do.call(tree_nj,
                           c(list(indel_distance(mat, missing = missing)), extra)),
      parsimony  = do.call(tree_parsimony, c(list(mat), extra)),
      iqtree     = do.call(tree_iqtree, c(list(mat), extra)),
      mix        = do.call(tree_mix, c(list(mat), extra)),
      vine       = do.call(tree_vine, c(list(mat), extra)),
      cassiopeia = do.call(tree_cassiopeia, c(list(mat), extra)))
  }

  trees <- stats::setNames(vector("list", length(methods)), methods)
  for (m in methods) {
    tree <- tryCatch(build_one(m), error = function(e) {
      if (identical(on_error, "stop")) {
        rlang::abort(
          sprintf("Tree method '%s' failed: %s", m, conditionMessage(e)),
          parent = e)
      }
      rlang::warn(sprintf("Skipping tree method '%s': %s", m, conditionMessage(e)))
      NULL
    })
    trees[[m]] <- tree
    if (!is.null(tree) && !is.null(output_dir)) {
      .require_tree_package("ape", "build_trees")
      ape::write.tree(tree, file = file.path(output_dir, sprintf("tree_%s.nwk", m)))
    }
  }

  ok <- names(trees)[!vapply(trees, is.null, logical(1L))]
  rlang::inform(sprintf(
    "build_trees: %d of %d method(s) succeeded%s.",
    length(ok), length(methods),
    if (length(ok)) paste0(" (", paste(ok, collapse = ", "), ")") else ""))
  trees
}
