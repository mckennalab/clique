#' Cassiopeia tree reconstruction (greedy / ILP / hybrid)
#'
#' Wraps the Python [Cassiopeia](https://github.com/YosefLab/Cassiopeia) package
#' via `reticulate`. Cassiopeia is an optional dependency: this function errors
#' with installation instructions if it isn't available.
#'
#' @param indel_mat Output `$matrix` from [build_indel_matrix()].
#' @param solver One of `"greedy"`, `"ilp"`, `"hybrid"`, `"upgma"`,
#'   `"neighbor_joining"`.
#' @param priors Optional named list of per-site indel priors.
#' @param missing Integer value used to encode missing data in `indel_mat`.
#' @param ... Named keyword arguments passed to the chosen Cassiopeia solver's
#'   constructor. For `solver = "hybrid"`, `top_solver` and `bottom_solver` may
#'   be supplied as Python solver objects; the defaults are VanillaGreedySolver
#'   and ILPSolver, with `cell_cutoff = 10`.
#' @return An `ape::phylo` object with original cell names as tip labels.
#' @examples
#' \dontrun{
#' tree <- tree_cassiopeia(im$matrix, solver = "greedy")
#' tree <- tree_cassiopeia(im$matrix, solver = "hybrid", priors = my_priors)
#' }
#' @export
tree_cassiopeia <- function(indel_mat,
                            solver = c("greedy", "ilp", "hybrid",
                                       "upgma", "neighbor_joining"),
                            priors = NULL,
                            missing = -1L,
                            ...) {
  solver <- match.arg(solver)
  mat <- .validate_indel_matrix(indel_mat, missing = missing)
  .require_tree_package("ape", "tree_cassiopeia")
  .require_tree_package(
    "reticulate", "tree_cassiopeia",
    "install.packages('reticulate'), then install Cassiopeia in its Python environment."
  )
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    rlang::abort("All arguments in `...` must be named Cassiopeia constructor arguments.")
  }

  cass <- tryCatch(
    suppressWarnings(reticulate::import(
      "cassiopeia", delay_load = FALSE, convert = FALSE
    )),
    error = function(e) {
      rlang::abort(c(
        "Python package 'cassiopeia' is required for tree_cassiopeia().",
        i = "Install it in the Python environment used by reticulate, for example with `pip install cassiopeia-lineage`.",
        i = conditionMessage(e)
      ))
    }
  )

  tip_map <- .tree_tip_map(rownames(mat))
  frame <- as.data.frame(mat, stringsAsFactors = FALSE)
  rownames(frame) <- names(tip_map)
  names(frame) <- as.character(seq_len(ncol(frame)) - 1L)
  py_matrix <- reticulate::r_to_py(frame, convert = TRUE)
  py_priors <- .cassiopeia_priors(priors, mat)
  tree_args <- list(
    character_matrix = py_matrix,
    missing_state_indicator = as.integer(missing)
  )
  if (!is.null(py_priors)) tree_args$priors <- py_priors
  cass_tree <- do.call(cass$data$CassiopeiaTree, tree_args)

  solver_object <- switch(
    solver,
    greedy = do.call(cass$solver$VanillaGreedySolver, dots),
    ilp = do.call(cass$solver$ILPSolver, dots),
    upgma = do.call(cass$solver$UPGMASolver, dots),
    neighbor_joining = do.call(cass$solver$NeighborJoiningSolver, dots),
    hybrid = {
      if (is.null(dots$top_solver)) {
        dots$top_solver <- cass$solver$VanillaGreedySolver()
      }
      if (is.null(dots$bottom_solver)) {
        dots$bottom_solver <- cass$solver$ILPSolver()
      }
      if (is.null(dots$cell_cutoff) && is.null(dots$lca_cutoff)) {
        dots$cell_cutoff <- 10L
      }
      do.call(cass$solver$HybridSolver, dots)
    }
  )

  tryCatch(
    solver_object$solve(cass_tree),
    error = function(e) {
      rlang::abort(c(
        sprintf("Cassiopeia's %s solver failed.", solver),
        i = conditionMessage(e),
        i = if (solver %in% c("ilp", "hybrid")) {
          "The ILP backend also requires a working Gurobi installation and license."
        } else NULL
      ))
    }
  )
  newick <- reticulate::py_to_r(cass_tree$get_newick(record_branch_lengths = TRUE))
  tree <- .read_backend_tree(text = newick, backend = "Cassiopeia")
  tree <- .restore_tree_tips(tree, tip_map, "Cassiopeia")
  .annotate_backend_tree(tree, paste0("cassiopeia-", solver))
}

#' @noRd
.cassiopeia_priors <- function(priors, mat) {
  if (is.null(priors)) return(NULL)
  if (!is.list(priors)) {
    rlang::abort("`priors` must be NULL or a list with one element per site.")
  }
  if (!is.null(names(priors))) {
    idx <- match(colnames(mat), names(priors))
    if (anyNA(idx)) {
      rlang::abort("Named `priors` must contain every `indel_mat` column name.")
    }
    priors <- priors[idx]
  } else if (length(priors) != ncol(mat)) {
    rlang::abort(sprintf("`priors` must have one element for each of %d sites.", ncol(mat)))
  }

  outer <- reticulate::dict(convert = FALSE)
  for (j in seq_along(priors)) {
    site <- priors[[j]]
    if (is.null(site)) next
    if (!is.numeric(site) || is.null(names(site)) || anyNA(site) ||
        any(!is.finite(site)) || any(site <= 0) || any(site > 1)) {
      rlang::abort(sprintf(
        "Prior for site '%s' must be a named numeric vector of probabilities in (0, 1].",
        colnames(mat)[j]
      ))
    }
    states <- suppressWarnings(as.integer(names(site)))
    if (anyNA(states) || any(as.character(states) != names(site)) || any(states < 0L)) {
      rlang::abort(sprintf("Prior names for site '%s' must be non-negative integer states.", colnames(mat)[j]))
    }
    inner <- reticulate::dict(convert = FALSE)
    for (k in seq_along(site)) {
      reticulate::py_set_item(inner, states[k], as.numeric(site[k]))
    }
    reticulate::py_set_item(outer, as.integer(j - 1L), inner)
  }
  outer
}
