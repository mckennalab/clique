#' Bayesian lineage-tree reconstruction with VINE
#'
#' Writes clique lineage characters in VINE's native CRISPR matrix format,
#' runs variational inference, and returns the posterior mean tree. VINE models
#' state `0` as unedited, positive integers as edited states, and `-1` as
#' missing. Its CRISPR mode is ultrametric.
#'
#' @details
#' Set `nj_only = TRUE` for VINE's fast initialization tree without variational
#' inference. Otherwise, VINE writes posterior tree samples to standard output
#' and its `--mean` tree is returned. With `return_posterior = TRUE`, parsed
#' samples are attached as the `clique_posterior` attribute.
#'
#' VINE is an optional external dependency available from Bioconda or Homebrew.
#'
#' @param indel_mat Output `$matrix` from [build_indel_matrix()].
#' @param model_type CRISPR rate-matrix parameterization: one matrix per site or
#'   one global matrix.
#' @param mutation_prior Uniform or empirical mutation-state frequencies.
#' @param dimensionality Embedding dimension. `NULL` uses VINE's logarithmic
#'   default, capped at one less than the number of cells so small data sets
#'   remain valid.
#' @param threads Positive thread count.
#' @param nsamples Number of posterior trees emitted after convergence.
#' @param nj_only If `TRUE`, return VINE's initial NJ/UPGMA tree without
#'   variational inference.
#' @param return_posterior Attach parsed posterior samples to the returned mean
#'   tree. Ignored for `nj_only = TRUE`.
#' @param binary Optional VINE executable name or path.
#' @param output_prefix Optional prefix for retained matrix, log, mean-tree, and
#'   posterior-sample files. Temporary files are removed when this is `NULL`.
#' @param extra_args Additional VINE arguments not managed by the wrapper.
#' @param missing Integer missing-data value. VINE requires `-1`.
#' @param timeout Process timeout in seconds (`Inf` for none).
#' @return An `ape::phylo` with original cell names as tip labels.
#' @examples
#' \dontrun{
#' tree <- tree_vine(im$matrix, threads = 8, nsamples = 500)
#' fast_tree <- tree_vine(im$matrix, nj_only = TRUE)
#' }
#' @export
tree_vine <- function(indel_mat,
                      model_type = c("sitewise", "global"),
                      mutation_prior = c("uniform", "empirical"),
                      dimensionality = NULL,
                      threads = 1L,
                      nsamples = 100L,
                      nj_only = FALSE,
                      return_posterior = FALSE,
                      binary = NULL,
                      output_prefix = NULL,
                      extra_args = character(),
                      missing = -1L,
                      timeout = Inf) {
  model_type <- match.arg(model_type)
  mutation_prior <- match.arg(mutation_prior)
  if (!identical(as.integer(missing), -1L)) {
    rlang::abort("VINE's CRISPR input format requires `missing = -1`.")
  }
  mat <- .validate_indel_matrix(indel_mat, missing = missing)
  if (is.null(dimensionality)) {
    dimensionality <- min(
      as.integer(round(3.25 + 0.92 * log(nrow(mat)))),
      nrow(mat) - 1L
    )
  }
  if (length(dimensionality) != 1L || is.na(dimensionality) ||
      dimensionality != as.integer(dimensionality) || dimensionality < 1L ||
      dimensionality >= nrow(mat)) {
    rlang::abort("`dimensionality` must be a positive integer smaller than the number of cells.")
  }
  if (length(threads) != 1L || is.na(threads) || threads != as.integer(threads) || threads < 1L) {
    rlang::abort("`threads` must be one positive integer.")
  }
  if (length(nsamples) != 1L || is.na(nsamples) || nsamples != as.integer(nsamples) || nsamples < 1L) {
    rlang::abort("`nsamples` must be one positive integer.")
  }
  if (!is.logical(nj_only) || length(nj_only) != 1L || is.na(nj_only) ||
      !is.logical(return_posterior) || length(return_posterior) != 1L || is.na(return_posterior)) {
    rlang::abort("`nj_only` and `return_posterior` must each be TRUE or FALSE.")
  }
  extra_args <- .validate_external_args(
    extra_args,
    c("--format", "-i", "--parallel", "-j", "--mean", "-m", "--nsamples", "-s",
      "--nj-only", "-0", "--crispr-modtype", "-Y", "--crispr-mutprior", "-p",
      "--logfile", "-l", "--dimensionality", "-D")
  )
  executable <- .resolve_tree_executable(
    binary, "vine", "VINE",
    "Install VINE with `conda install -c conda-forge -c bioconda vine-phylo` or `brew install CshlSiepelLab/tools/vine`."
  )

  tip_map <- .tree_tip_map(rownames(mat))
  work <- .tree_workspace(output_prefix, "vine")
  if (work$temporary) on.exit(unlink(work$directory, recursive = TRUE), add = TRUE)
  input_file <- paste0(work$prefix, ".vine.tsv")
  mean_file <- paste0(work$prefix, ".vine.mean.nwk")
  log_file <- paste0(work$prefix, ".vine.log")
  samples_file <- paste0(work$prefix, ".vine.samples.nwk")
  .write_vine_matrix(mat, input_file, tip_map)

  args <- c(
    input_file,
    "--format", "CRISPR",
    "--dimensionality", as.character(as.integer(dimensionality)),
    "--parallel", as.character(as.integer(threads)),
    "--crispr-modtype", toupper(model_type),
    "--crispr-mutprior", if (mutation_prior == "uniform") "UNIF" else "EMPIRICAL",
    "--silent"
  )
  if (nj_only) {
    args <- c(args, "--nj-only")
  } else {
    args <- c(
      args,
      "--mean", mean_file,
      "--nsamples", as.character(as.integer(nsamples)),
      "--logfile", log_file
    )
  }
  args <- c(args, extra_args)
  result <- .run_tree_process(
    executable, args, "VINE", timeout = timeout, abort_on_error = FALSE
  )

  newick_lines <- .vine_newick_lines(result$stdout)
  if (!work$temporary) writeLines(newick_lines, samples_file, useBytes = TRUE)
  if (nj_only) {
    if (length(newick_lines) == 0L) {
      .abort_tree_process(result, executable, args, "VINE")
    }
    tree <- .read_backend_tree(text = newick_lines[1L], backend = "VINE")
  } else {
    if (!file.exists(mean_file) || is.na(file.info(mean_file)$size) || file.info(mean_file)$size == 0L) {
      .abort_tree_process(result, executable, args, "VINE")
    }
    tree <- .read_backend_tree(path = mean_file, backend = "VINE")
  }
  if (!identical(result$status, 0L)) {
    rlang::warn(sprintf(
      "VINE returned a usable tree but exited with status %s during shutdown.",
      result$status
    ))
  }
  tree <- .restore_tree_tips(tree, tip_map, "VINE")

  if (!nj_only && return_posterior) {
    if (length(newick_lines) == 0L) {
      rlang::abort("VINE completed without writing posterior tree samples.")
    }
    posterior <- .read_backend_tree(
      text = paste(newick_lines, collapse = "\n"), backend = "VINE", all = TRUE
    )
    posterior <- .restore_tree_tips(posterior, tip_map, "VINE")
    attr(tree, "clique_posterior") <- posterior
  }
  files <- if (!work$temporary) {
    out <- c(matrix = input_file, samples = samples_file)
    if (!nj_only) out <- c(out, mean = mean_file, log = log_file)
    out
  } else NULL
  .annotate_backend_tree(tree, "vine", c(executable, args), files)
}

#' @noRd
.write_vine_matrix <- function(mat, path, tip_map) {
  out <- data.frame(cell = names(tip_map), mat, check.names = FALSE)
  names(out)[-1L] <- colnames(mat)
  utils::write.table(
    out, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "-1"
  )
  invisible(path)
}

#' @noRd
.vine_newick_lines <- function(stdout) {
  if (is.null(stdout) || !nzchar(stdout)) return(character())
  lines <- trimws(strsplit(stdout, "\n", fixed = TRUE)[[1L]])
  lines[nzchar(lines) & grepl("(", lines, fixed = TRUE) & grepl(";$", lines)]
}
