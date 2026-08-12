# Shared validation and process helpers for lineage-tree backends.

#' @noRd
.validate_indel_matrix <- function(indel_mat, missing = -1L, min_cells = 3L) {
  if (inherits(indel_mat, "Matrix")) indel_mat <- as.matrix(indel_mat)
  if (is.data.frame(indel_mat)) indel_mat <- as.matrix(indel_mat)
  if (!is.matrix(indel_mat) || !is.numeric(indel_mat)) {
    rlang::abort("`indel_mat` must be a numeric matrix with cells in rows and sites in columns.")
  }
  if (length(missing) != 1L || is.na(missing) || missing != as.integer(missing)) {
    rlang::abort("`missing` must be one integer value.")
  }
  missing <- as.integer(missing)
  if (nrow(indel_mat) < min_cells) {
    rlang::abort(sprintf("Tree reconstruction needs at least %d cells.", min_cells))
  }
  if (ncol(indel_mat) < 1L) {
    rlang::abort("`indel_mat` must contain at least one site.")
  }
  if (anyNA(indel_mat) || any(!is.finite(indel_mat))) {
    rlang::abort(sprintf(
      "`indel_mat` contains NA or non-finite values; encode missing calls as %d.",
      missing
    ))
  }
  if (any(indel_mat != round(indel_mat))) {
    rlang::abort("All character states in `indel_mat` must be integers.")
  }

  indel_mat <- matrix(
    as.integer(indel_mat),
    nrow = nrow(indel_mat),
    ncol = ncol(indel_mat),
    dimnames = dimnames(indel_mat)
  )
  bad <- indel_mat < 0L & indel_mat != missing
  if (any(bad)) {
    rlang::abort(sprintf(
      "Character states must be non-negative integers or the missing value %d.",
      missing
    ))
  }

  if (is.null(rownames(indel_mat))) {
    rownames(indel_mat) <- paste0("cell", seq_len(nrow(indel_mat)))
  }
  if (anyNA(rownames(indel_mat)) || any(!nzchar(rownames(indel_mat)))) {
    rlang::abort("Cell names in `rownames(indel_mat)` must be non-empty and non-missing.")
  }
  if (anyDuplicated(rownames(indel_mat))) {
    rlang::abort("Cell names in `rownames(indel_mat)` must be unique.")
  }
  if (is.null(colnames(indel_mat))) {
    colnames(indel_mat) <- paste0("site", seq_len(ncol(indel_mat)))
  }

  all_missing <- colSums(indel_mat != missing) == 0L
  if (all(all_missing)) {
    rlang::abort("Every site in `indel_mat` is entirely missing.")
  }
  if (any(all_missing)) {
    rlang::warn(sprintf(
      "Dropping %d site(s) that are missing in every cell.",
      sum(all_missing)
    ))
    indel_mat <- indel_mat[, !all_missing, drop = FALSE]
  }
  indel_mat
}

#' @noRd
.tree_tip_map <- function(labels) {
  if (length(labels) > 9999999L) {
    rlang::abort("External tree backends currently support at most 9,999,999 cells.")
  }
  safe <- sprintf("clq%07d", seq_along(labels))
  stats::setNames(as.character(labels), safe)
}

#' @noRd
.restore_tree_tips <- function(tree, tip_map, backend) {
  restore_one <- function(x) {
    idx <- match(x$tip.label, names(tip_map))
    if (anyNA(idx)) {
      unexpected <- unique(x$tip.label[is.na(idx)])
      rlang::abort(c(
        sprintf("%s returned unexpected tip labels.", backend),
        i = paste(utils::head(unexpected, 5L), collapse = ", ")
      ))
    }
    if (length(idx) != length(tip_map) || anyDuplicated(idx)) {
      rlang::abort(sprintf(
        "%s returned %d unique tips; expected %d.",
        backend, length(unique(idx)), length(tip_map)
      ))
    }
    x$tip.label <- unname(tip_map[idx])
    x
  }

  if (inherits(tree, "multiPhylo")) {
    tree[] <- lapply(unclass(tree), restore_one)
    return(tree)
  }
  restore_one(tree)
}

#' @noRd
.require_tree_package <- function(package, caller, install = NULL) {
  if (requireNamespace(package, quietly = TRUE)) return(invisible(TRUE))
  if (is.null(install)) install <- sprintf("install.packages('%s')", package)
  rlang::abort(c(
    sprintf("Package '%s' is required for %s().", package, caller),
    i = install
  ))
}

#' @noRd
.resolve_tree_executable <- function(binary, candidates, backend, install) {
  if (!is.null(binary)) {
    if (!is.character(binary) || length(binary) != 1L || !nzchar(binary)) {
      rlang::abort("`binary` must be NULL or one non-empty executable name/path.")
    }
    found <- Sys.which(binary)
    if (!nzchar(found) && file.exists(binary) && file.access(binary, 1L) == 0L) {
      found <- normalizePath(binary, mustWork = TRUE)
    }
  } else {
    found <- Sys.which(candidates)
    found <- unname(found[nzchar(found)][1L])
  }

  if (length(found) == 0L || is.na(found) || !nzchar(found)) {
    rlang::abort(c(
      sprintf("Could not find the %s executable.", backend),
      i = install,
      i = "Pass its executable name or full path with `binary`."
    ))
  }
  normalizePath(found, mustWork = TRUE)
}

#' @noRd
.validate_external_args <- function(extra_args, reserved = character()) {
  if (!is.character(extra_args) || anyNA(extra_args)) {
    rlang::abort("`extra_args` must be a character vector without NA values.")
  }
  conflict <- intersect(tolower(extra_args), tolower(reserved))
  if (length(conflict)) {
    rlang::abort(sprintf(
      "`extra_args` cannot override wrapper-managed option(s): %s.",
      paste(unique(conflict), collapse = ", ")
    ))
  }
  extra_args
}

#' @noRd
.run_tree_process <- function(command,
                              args = character(),
                              backend,
                              wd = NULL,
                              timeout = Inf,
                              stdin = NULL,
                              abort_on_error = TRUE) {
  result <- processx::run(
    command = command,
    args = args,
    wd = wd,
    echo = FALSE,
    echo_cmd = FALSE,
    spinner = FALSE,
    error_on_status = FALSE,
    timeout = timeout,
    stdin = stdin
  )
  if (!identical(result$status, 0L) && abort_on_error) {
    .abort_tree_process(result, command, args, backend)
  }
  result
}

#' @noRd
.abort_tree_process <- function(result, command, args, backend) {
  stderr <- if (is.null(result$stderr)) "" else result$stderr
  details <- paste(utils::tail(strsplit(stderr, "\n", fixed = TRUE)[[1L]], 20L),
                   collapse = "\n")
  rlang::abort(c(
    sprintf("%s failed with exit status %s.", backend, result$status),
    i = paste(c(command, args), collapse = " "),
    i = details
  ))
}

#' @noRd
.tree_workspace <- function(output_prefix, backend) {
  if (is.null(output_prefix)) {
    directory <- tempfile(paste0("cliqueR-", backend, "-"))
    dir.create(directory, recursive = TRUE)
    return(list(
      directory = directory,
      prefix = file.path(directory, "run"),
      temporary = TRUE
    ))
  }
  if (!is.character(output_prefix) || length(output_prefix) != 1L ||
      is.na(output_prefix) || !nzchar(output_prefix)) {
    rlang::abort("`output_prefix` must be NULL or one non-empty path.")
  }
  prefix <- normalizePath(path.expand(output_prefix), mustWork = FALSE)
  directory <- dirname(prefix)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
    rlang::abort(sprintf("Could not create output directory '%s'.", directory))
  }
  list(directory = directory, prefix = prefix, temporary = FALSE)
}

#' @noRd
.read_backend_tree <- function(path = NULL, text = NULL, backend, all = FALSE) {
  .require_tree_package("ape", backend)
  tree <- tryCatch(
    {
      if (!is.null(path)) ape::read.tree(path) else ape::read.tree(text = text)
    },
    error = function(e) {
      rlang::abort(c(
        sprintf("Could not parse the Newick tree returned by %s.", backend),
        i = conditionMessage(e)
      ))
    }
  )
  if (is.null(tree) || length(tree) == 0L) {
    rlang::abort(sprintf("%s did not return a Newick tree.", backend))
  }
  if (!all && inherits(tree, "multiPhylo")) tree <- tree[[1L]]
  tree
}

#' @noRd
.annotate_backend_tree <- function(tree, backend, command = NULL, files = NULL) {
  attr(tree, "clique_backend") <- backend
  if (!is.null(command)) attr(tree, "clique_command") <- command
  if (!is.null(files)) attr(tree, "clique_files") <- files
  tree
}
