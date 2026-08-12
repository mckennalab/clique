#!/usr/bin/env Rscript

# Introductory cliqueR tree-building workflow.
#
# Run with the bundled example:
#   Rscript inst/examples/tree-building-introduction.R
#
# Run with your own reference tree and clique character matrix:
#   Rscript inst/examples/tree-building-introduction.R \
#     --tree=reference.nwk \
#     --matrix=characters.tsv \
#     --methods=nj,parsimony,iqtree \
#     --output=tree-building-results

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript tree-building-introduction.R [options]\n\n",
    "Options:\n",
    "  --tree=PATH       Reference Newick tree.\n",
    "  --matrix=PATH     Tab-separated character matrix; first column is the cell ID.\n",
    "  --methods=LIST    Comma-separated cliqueR methods (default: nj,parsimony).\n",
    "  --output=DIR      Output directory (default: tree-building-results).\n",
    "  --help            Show this help.\n\n",
    "Character states must be integers: 0 is unedited, positive values are edit\n",
    "alleles, and -1 is missing. Tree tips and matrix cell IDs must match.\n"
  ))
  quit(save = "no", status = status)
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0L) return(normalizePath(getwd()))
  normalizePath(sub("^--file=", "", file_arg[[1L]]))
}

parse_options <- function(args, example_dir) {
  options <- list(
    tree = file.path(example_dir, "data", "tree-building-reference.nwk"),
    matrix = file.path(example_dir, "data", "tree-building-characters.tsv"),
    methods = c("nj", "parsimony"),
    output = file.path(getwd(), "tree-building-results")
  )
  if (any(args %in% c("-h", "--help"))) usage()

  for (arg in args) {
    if (!grepl("^--[^=]+=.+$", arg)) {
      stop("Unknown argument '", arg, "'. Use --help for usage.", call. = FALSE)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!key %in% names(options)) {
      stop("Unknown option '--", key, "'. Use --help for usage.", call. = FALSE)
    }
    options[[key]] <- value
  }

  options$methods <- trimws(unlist(
    strsplit(options$methods, ",", fixed = TRUE),
    use.names = FALSE
  ))
  options$methods <- options$methods[nzchar(options$methods)]
  if (length(options$methods) == 0L) {
    stop("--methods must contain at least one method.", call. = FALSE)
  }
  options
}

read_character_matrix <- function(path) {
  if (!file.exists(path)) stop("Character matrix not found: ", path, call. = FALSE)
  input <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  if (ncol(input) < 2L) {
    stop("The character matrix needs a cell-ID column and at least one site.", call. = FALSE)
  }

  cell_ids <- as.character(input[[1L]])
  if (anyNA(cell_ids) || any(!nzchar(cell_ids)) || anyDuplicated(cell_ids)) {
    stop("The first matrix column must contain unique, non-empty cell IDs.", call. = FALSE)
  }

  values <- as.matrix(input[-1L])
  numeric_values <- suppressWarnings(matrix(
    as.numeric(values),
    nrow = nrow(values),
    ncol = ncol(values),
    dimnames = list(cell_ids, names(input)[-1L])
  ))
  invalid <- is.na(numeric_values) & !is.na(values)
  if (any(invalid)) {
    stop("Every character state must be an integer, NA, or blank.", call. = FALSE)
  }
  if (any(numeric_values != round(numeric_values), na.rm = TRUE)) {
    stop("Every observed character state must be an integer.", call. = FALSE)
  }
  numeric_values[is.na(numeric_values)] <- -1
  storage.mode(numeric_values) <- "integer"
  numeric_values
}

first_tree <- function(tree) {
  if (inherits(tree, "multiPhylo")) tree[[1L]] else tree
}

character_data <- function(indel_mat) {
  states <- sort(unique(indel_mat[indel_mat != -1L]))
  encoded <- matrix(
    as.character(indel_mat),
    nrow(indel_mat),
    ncol(indel_mat),
    dimnames = dimnames(indel_mat)
  )
  encoded[indel_mat == -1L] <- "?"
  phangorn::phyDat(encoded, type = "USER", levels = as.character(states))
}

unit_edge_lengths <- function(tree) {
  if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
  tree
}

tip_distance_correlation <- function(reference, candidate) {
  reference <- unit_edge_lengths(reference)
  candidate <- unit_edge_lengths(candidate)
  reference_dist <- ape::cophenetic.phylo(reference)
  candidate_dist <- ape::cophenetic.phylo(candidate)[rownames(reference_dist), colnames(reference_dist)]
  keep <- upper.tri(reference_dist)
  suppressWarnings(stats::cor(
    reference_dist[keep], candidate_dist[keep], method = "spearman"
  ))
}

compare_to_reference <- function(reference, trees, indel_mat) {
  data <- character_data(indel_mat)
  candidates <- c(list(reference = reference), trees)
  candidates <- lapply(candidates, first_tree)

  rows <- lapply(names(candidates), function(method) {
    tree <- candidates[[method]]
    data.frame(
      method = method,
      backend = if (method == "reference") "input-newick" else
        or_else(attr(tree, "clique_backend"), method),
      tips = length(tree$tip.label),
      rf_distance = as.numeric(phangorn::RF.dist(
        reference, tree, normalize = FALSE, rooted = FALSE
      )),
      normalized_rf = as.numeric(phangorn::RF.dist(
        reference, tree, normalize = TRUE, rooted = FALSE
      )),
      parsimony_score = as.numeric(phangorn::parsimony(tree, data)),
      tip_distance_spearman = tip_distance_correlation(reference, tree),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pairwise_rf <- function(reference, trees) {
  candidates <- c(list(reference = reference), trees)
  candidates <- lapply(candidates, first_tree)
  labels <- names(candidates)
  result <- matrix(0, length(candidates), length(candidates), dimnames = list(labels, labels))
  for (i in seq_along(candidates)) {
    if (i == length(candidates)) break
    for (j in (i + 1L):length(candidates)) {
      distance <- as.numeric(phangorn::RF.dist(
        candidates[[i]], candidates[[j]], normalize = TRUE, rooted = FALSE
      ))
      result[i, j] <- result[j, i] <- distance
    }
  }
  result
}

plot_trees <- function(reference, trees, path) {
  candidates <- c(list(reference = reference), trees)
  candidates <- lapply(candidates, first_tree)
  columns <- min(2L, length(candidates))
  rows <- ceiling(length(candidates) / columns)
  grDevices::pdf(path, width = 6 * columns, height = 4.5 * rows)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(rows, columns), mar = c(1, 1, 3, 1))
  for (method in names(candidates)) {
    ape::plot.phylo(candidates[[method]], cex = 0.8, no.margin = TRUE)
    graphics::title(main = method)
  }
}

or_else <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

main <- function() {
  if (!requireNamespace("cliqueR", quietly = TRUE)) {
    stop("Install cliqueR before running this example.", call. = FALSE)
  }
  for (package in c("ape", "phangorn")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the suggested R package '", package, "' to run this example.", call. = FALSE)
    }
  }

  example_dir <- dirname(script_path())
  options <- parse_options(commandArgs(trailingOnly = TRUE), example_dir)
  if (!file.exists(options$tree)) stop("Reference tree not found: ", options$tree, call. = FALSE)

  reference <- ape::read.tree(options$tree)
  if (inherits(reference, "multiPhylo") || !inherits(reference, "phylo")) {
    stop("--tree must contain exactly one Newick tree.", call. = FALSE)
  }
  if (anyDuplicated(reference$tip.label)) {
    stop("Reference-tree tip labels must be unique.", call. = FALSE)
  }

  indel_mat <- read_character_matrix(options$matrix)
  missing_tips <- setdiff(reference$tip.label, rownames(indel_mat))
  extra_cells <- setdiff(rownames(indel_mat), reference$tip.label)
  if (length(missing_tips) || length(extra_cells)) {
    stop(
      "Reference tips and matrix cells differ. Missing from matrix: ",
      paste(missing_tips, collapse = ", "),
      "; absent from tree: ", paste(extra_cells, collapse = ", "),
      call. = FALSE
    )
  }
  indel_mat <- indel_mat[reference$tip.label, , drop = FALSE]

  dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
  ape::write.tree(reference, file.path(options$output, "tree_reference.nwk"))
  utils::write.table(
    data.frame(cell_id = rownames(indel_mat), indel_mat, check.names = FALSE),
    file.path(options$output, "character_matrix.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  backend_args <- list(
    iqtree = list(seed = 1L),
    vine = list(nj_only = TRUE),
    cassiopeia = list(solver = "greedy")
  )
  trees <- cliqueR::build_trees(
    indel_mat,
    methods = options$methods,
    output_dir = options$output,
    args = backend_args,
    on_error = "warn"
  )
  successful <- trees[!vapply(trees, is.null, logical(1L))]
  if (length(successful) == 0L) {
    stop("No tree-building method succeeded.", call. = FALSE)
  }

  comparison <- compare_to_reference(reference, successful, indel_mat)
  pairwise <- pairwise_rf(reference, successful)
  utils::write.table(
    comparison,
    file.path(options$output, "tree_comparison.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    data.frame(tree = rownames(pairwise), pairwise, check.names = FALSE),
    file.path(options$output, "pairwise_normalized_rf.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  plot_trees(reference, successful, file.path(options$output, "tree_comparison.pdf"))

  cat("\nTree comparison (lower RF and parsimony; higher correlation are better):\n")
  print(comparison, row.names = FALSE, digits = 3)
  failed <- setdiff(options$methods, names(successful))
  if (length(failed)) cat("\nSkipped methods: ", paste(failed, collapse = ", "), "\n", sep = "")
  cat("\nResults written to ", normalizePath(options$output), "\n", sep = "")
}

if (sys.nframe() == 0L) main()
