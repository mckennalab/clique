#' Maximum-likelihood lineage tree with IQ-TREE 2
#'
#' Encodes a clique integer character matrix as a NEXUS morphological
#' alignment, runs IQ-TREE 2, and imports the best tree as an `ape::phylo`.
#' State labels are recoded independently at each site because edit-state
#' integers are categorical, not ordered measurements.
#'
#' @details
#' IQ-TREE supports at most 32 morphological states (`0`-`9`, `A`-`V`) at one
#' site. This function stops rather than merging alleles when a site exceeds
#' that limit. The default `MK+FQ` model does not use ascertainment correction
#' because clique matrices may contain invariant sites. If invariant sites were
#' removed before calling this function, consider `model = "MK+FQ+ASC"`.
#'
#' IQ-TREE is an optional external dependency. The wrapper searches for
#' `iqtree2`, then `iqtree`, unless `binary` is supplied.
#'
#' @param indel_mat Output `$matrix` from [build_indel_matrix()].
#' @param model IQ-TREE morphological model.
#' @param threads Positive thread count or `"AUTO"`.
#' @param seed Optional IQ-TREE random seed.
#' @param outgroup Optional cell name from `rownames(indel_mat)`.
#' @param binary Optional IQ-TREE executable name or path.
#' @param output_prefix Optional prefix for retained NEXUS and IQ-TREE output
#'   files. Temporary files are removed when this is `NULL`.
#' @param extra_args Additional IQ-TREE command-line arguments. Arguments that
#'   replace wrapper-managed input, output, model, thread, seed, or outgroup
#'   options are rejected.
#' @param missing Integer value used for missing data.
#' @param timeout Process timeout in seconds (`Inf` for none).
#' @return An `ape::phylo` with original cell names as tip labels. Attributes
#'   `clique_backend`, `clique_command`, and, when retained, `clique_files`
#'   record execution details.
#' @examples
#' \dontrun{
#' tree <- tree_iqtree(im$matrix, threads = 8, model = "MK+FQ")
#' }
#' @export
tree_iqtree <- function(indel_mat,
                        model = "MK+FQ",
                        threads = 1L,
                        seed = 1L,
                        outgroup = NULL,
                        binary = NULL,
                        output_prefix = NULL,
                        extra_args = character(),
                        missing = -1L,
                        timeout = Inf) {
  mat <- .validate_indel_matrix(indel_mat, missing = missing)
  if (!is.character(model) || length(model) != 1L || !nzchar(model)) {
    rlang::abort("`model` must be one non-empty IQ-TREE model string.")
  }
  if (length(threads) != 1L || is.na(threads) ||
      !(identical(toupper(as.character(threads)), "AUTO") ||
        (!is.na(suppressWarnings(as.integer(threads))) &&
         as.numeric(threads) == as.integer(threads) && as.integer(threads) >= 1L))) {
    rlang::abort("`threads` must be a positive integer or 'AUTO'.")
  }
  threads <- if (toupper(as.character(threads)) == "AUTO") "AUTO" else as.character(as.integer(threads))
  if (!is.null(seed) && (length(seed) != 1L || is.na(seed) || seed != as.integer(seed))) {
    rlang::abort("`seed` must be NULL or one integer.")
  }

  tip_map <- .tree_tip_map(rownames(mat))
  outgroup_safe <- NULL
  if (!is.null(outgroup)) {
    if (!is.character(outgroup) || length(outgroup) != 1L || !outgroup %in% unname(tip_map)) {
      rlang::abort("`outgroup` must name one cell in `rownames(indel_mat)`.")
    }
    outgroup_safe <- names(tip_map)[match(outgroup, unname(tip_map))]
  }
  extra_args <- .validate_external_args(
    extra_args,
    c("-s", "--seqfile", "-st", "-m", "-nt", "-seed", "-o", "-pre", "--prefix")
  )
  executable <- .resolve_tree_executable(
    binary, c("iqtree2", "iqtree"), "IQ-TREE",
    "Install IQ-TREE 2 and ensure `iqtree2` is on PATH."
  )

  work <- .tree_workspace(output_prefix, "iqtree")
  if (work$temporary) on.exit(unlink(work$directory, recursive = TRUE), add = TRUE)
  nexus <- paste0(work$prefix, ".nex")
  .write_iqtree_nexus(mat, nexus, tip_map, missing)

  args <- c(
    "-s", nexus,
    "-st", "MORPH",
    "-m", model,
    "-nt", threads,
    "-pre", work$prefix,
    "-quiet"
  )
  if (!is.null(seed)) args <- c(args, "-seed", as.character(as.integer(seed)))
  if (!is.null(outgroup_safe)) args <- c(args, "-o", outgroup_safe)
  args <- c(args, extra_args)

  .run_tree_process(executable, args, "IQ-TREE", timeout = timeout)
  tree_file <- paste0(work$prefix, ".treefile")
  if (!file.exists(tree_file)) {
    rlang::abort(sprintf("IQ-TREE completed without creating '%s'.", tree_file))
  }
  tree <- .read_backend_tree(path = tree_file, backend = "IQ-TREE")
  tree <- .restore_tree_tips(tree, tip_map, "IQ-TREE")
  files <- if (!work$temporary) c(alignment = nexus, tree = tree_file) else NULL
  .annotate_backend_tree(tree, "iqtree2", c(executable, args), files)
}

#' @noRd
.write_iqtree_nexus <- function(mat, path, tip_map, missing) {
  symbols <- c(as.character(0:9), LETTERS[1:22])
  encoded <- matrix("?", nrow(mat), ncol(mat))
  max_states <- 0L
  for (j in seq_len(ncol(mat))) {
    observed <- sort(unique(mat[mat[, j] != missing, j]))
    observed <- c(intersect(0L, observed), setdiff(observed, 0L))
    if (length(observed) > length(symbols)) {
      rlang::abort(sprintf(
        "Site '%s' has %d observed states; IQ-TREE supports at most 32 morphological states.",
        colnames(mat)[j], length(observed)
      ))
    }
    codes <- stats::setNames(symbols[seq_along(observed)], observed)
    present <- mat[, j] != missing
    encoded[present, j] <- unname(codes[as.character(mat[present, j])])
    max_states <- max(max_states, length(observed))
  }

  lines <- c(
    "#NEXUS",
    "BEGIN DATA;",
    sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", nrow(mat), ncol(mat)),
    sprintf(
      "  FORMAT DATATYPE=STANDARD SYMBOLS=\"%s\" MISSING=? GAP=-;",
      paste0(symbols[seq_len(max(1L, max_states))], collapse = "")
    ),
    "  MATRIX",
    sprintf("  %-10s %s", names(tip_map), apply(encoded, 1L, paste0, collapse = "")),
    "  ;",
    "END;"
  )
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}
