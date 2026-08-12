#' Parsimony lineage tree with PHYLIP MIX
#'
#' Runs the PHYLIP `mix` program on a binary expansion of the clique character
#' matrix. Each positive edit state becomes a presence/absence character at its
#' site; state `0` is absent for every edit character, and missing calls become
#' `?`. This avoids incorrectly treating arbitrary edit-state integers as an
#' ordered scale.
#'
#' @details
#' Camin-Sokal parsimony is the default because it models irreversible gains of
#' edit alleles. Wagner parsimony permits gains and losses. PHYLIP can emit
#' multiple equally parsimonious trees; by default the first is returned and
#' the number found is recorded in `clique_tied_trees`. Set `all = TRUE` to
#' return all trees as an `ape::multiPhylo`.
#'
#' PHYLIP MIX is an optional external dependency. It uses fixed working-file
#' names, so this wrapper always runs it in an isolated temporary directory.
#'
#' @param indel_mat Output `$matrix` from [build_indel_matrix()].
#' @param method Camin-Sokal (irreversible) or Wagner parsimony.
#' @param outgroup Optional cell name from `rownames(indel_mat)`.
#' @param all Return all equally parsimonious trees instead of the first.
#' @param binary Optional path or name for the PHYLIP `mix` executable.
#' @param output_prefix Optional prefix under which the generated MIX input,
#'   menu, report, tree, and screen output are retained.
#' @param missing Integer value used for missing data.
#' @param timeout Process timeout in seconds (`Inf` for none).
#' @return An `ape::phylo`, or an `ape::multiPhylo` when `all = TRUE`.
#' @examples
#' \dontrun{
#' tree <- tree_mix(im$matrix, method = "camin-sokal")
#' }
#' @export
tree_mix <- function(indel_mat,
                     method = c("camin-sokal", "wagner"),
                     outgroup = NULL,
                     all = FALSE,
                     binary = NULL,
                     output_prefix = NULL,
                     missing = -1L,
                     timeout = Inf) {
  method <- match.arg(method)
  mat <- .validate_indel_matrix(indel_mat, missing = missing)
  if (!is.logical(all) || length(all) != 1L || is.na(all)) {
    rlang::abort("`all` must be TRUE or FALSE.")
  }
  outgroup_index <- NULL
  if (!is.null(outgroup)) {
    if (!is.character(outgroup) || length(outgroup) != 1L || !outgroup %in% rownames(mat)) {
      rlang::abort("`outgroup` must name one cell in `rownames(indel_mat)`.")
    }
    outgroup_index <- match(outgroup, rownames(mat))
  }
  executable <- .resolve_tree_executable(
    binary, "mix", "PHYLIP MIX",
    "Install PHYLIP and ensure its `mix` executable is on PATH."
  )
  tip_map <- .tree_tip_map(rownames(mat))
  binary_mat <- .mix_binary_characters(mat, missing)

  work_dir <- tempfile("cliqueR-mix-")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  infile <- file.path(work_dir, "infile")
  menu_file <- file.path(work_dir, "menu")
  .write_mix_input(binary_mat, infile, names(tip_map))

  menu <- character()
  if (method == "camin-sokal") menu <- c(menu, "P")
  if (!is.null(outgroup_index)) menu <- c(menu, "O", as.character(outgroup_index))
  writeLines(c(menu, "Y"), menu_file, useBytes = TRUE)
  result <- .run_tree_process(
    executable, backend = "PHYLIP MIX", wd = work_dir,
    timeout = timeout, stdin = menu_file
  )

  candidates <- file.path(work_dir, c("outtree", "treefile"))
  tree_file <- candidates[file.exists(candidates)][1L]
  if (length(tree_file) == 0L || is.na(tree_file)) {
    rlang::abort("PHYLIP MIX completed without creating an outtree/treefile.")
  }
  trees <- .read_backend_tree(path = tree_file, backend = "PHYLIP MIX", all = TRUE)
  if (!inherits(trees, "multiPhylo")) {
    trees <- structure(list(trees), class = "multiPhylo")
  }
  trees <- .restore_tree_tips(trees, tip_map, "PHYLIP MIX")
  n_trees <- length(trees)
  tree <- if (all) trees else trees[[1L]]
  attr(tree, "clique_tied_trees") <- n_trees

  files <- NULL
  if (!is.null(output_prefix)) {
    prefix <- .tree_workspace(output_prefix, "mix")$prefix
    files <- c(
      input = paste0(prefix, ".mix.infile"),
      menu = paste0(prefix, ".mix.menu"),
      report = paste0(prefix, ".mix.outfile"),
      tree = paste0(prefix, ".mix.outtree"),
      screen = paste0(prefix, ".mix.screenout")
    )
    file.copy(infile, files[["input"]], overwrite = TRUE)
    file.copy(menu_file, files[["menu"]], overwrite = TRUE)
    outfile <- file.path(work_dir, "outfile")
    if (file.exists(outfile)) file.copy(outfile, files[["report"]], overwrite = TRUE)
    file.copy(tree_file, files[["tree"]], overwrite = TRUE)
    writeLines(result$stdout, files[["screen"]], useBytes = TRUE)
    files <- files[file.exists(files)]
  }
  .annotate_backend_tree(tree, "phylip-mix", executable, files)
}

#' @noRd
.mix_binary_characters <- function(mat, missing) {
  columns <- list()
  labels <- character()
  for (j in seq_len(ncol(mat))) {
    states <- sort(unique(mat[mat[, j] != missing, j]))
    edited <- states[states != 0L]
    for (state in edited) {
      columns[[length(columns) + 1L]] <- ifelse(
        mat[, j] == missing, "?", ifelse(mat[, j] == state, "1", "0")
      )
      labels <- c(labels, sprintf("%s=%d", colnames(mat)[j], state))
    }
  }
  if (length(columns) == 0L) {
    rlang::abort("PHYLIP MIX needs at least one positive edit state.")
  }
  out <- do.call(cbind, columns)
  dimnames(out) <- list(rownames(mat), labels)
  out
}

#' @noRd
.write_mix_input <- function(binary_mat, path, safe_labels) {
  patterns <- apply(binary_mat, 1L, paste0, collapse = "")
  lines <- c(
    sprintf("%5d %5d", nrow(binary_mat), ncol(binary_mat)),
    sprintf("%-10s%s", safe_labels, patterns)
  )
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}
