# Interactive companion to tree-building-introduction.R.
#
# From the cliqueR package source directory:
#   source("inst/examples/tree-building-introduction-interactive.R")
#   result <- run_tree_building_example()
#   result$comparison
#   result$pairwise_rf
#
# With your own data:
#   result <- run_tree_building_example(
#     reference_tree = "reference.nwk",
#     character_matrix = "characters.tsv",
#     methods = c("nj", "parsimony", "iqtree"),
#     output_dir = "tree-building-results"
#   )

.tree_building_interactive_dir <- local({
  source_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (is.null(source_file)) {
    frame_files <- vapply(
      sys.frames(),
      function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile,
      character(1L)
    )
    frame_files <- frame_files[!is.na(frame_files)]
    if (length(frame_files)) source_file <- frame_files[[length(frame_files)]]
  }
  if (is.null(source_file)) getwd() else dirname(normalizePath(source_file))
})

.tree_building_helpers <- new.env(parent = baseenv())
.tree_building_cli <- file.path(
  .tree_building_interactive_dir,
  "tree-building-introduction.R"
)
if (!file.exists(.tree_building_cli)) {
  stop("Could not find tree-building-introduction.R beside this script.", call. = FALSE)
}
sys.source(.tree_building_cli, envir = .tree_building_helpers)

#' Locate the bundled introductory tree-building data
#'
#' @return A named character vector containing the reference Newick and
#'   character-matrix paths.
tree_building_example_files <- function() {
  paths <- c(
    reference_tree = file.path(
      .tree_building_interactive_dir,
      "data",
      "tree-building-reference.nwk"
    ),
    character_matrix = file.path(
      .tree_building_interactive_dir,
      "data",
      "tree-building-characters.tsv"
    )
  )
  if (any(!file.exists(paths))) {
    stop("The bundled tree-building example data could not be found.", call. = FALSE)
  }
  paths
}

.interactive_character_matrix <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    return(.tree_building_helpers$read_character_matrix(x))
  }

  if (is.data.frame(x)) {
    default_rows <- identical(rownames(x), as.character(seq_len(nrow(x))))
    if (default_rows && ncol(x) >= 2L &&
        (is.character(x[[1L]]) || is.factor(x[[1L]]))) {
      cell_ids <- as.character(x[[1L]])
      x <- x[-1L]
      rownames(x) <- cell_ids
    }
  }
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("character_matrix must be a TSV path, matrix, or data frame.", call. = FALSE)
  }

  values <- as.matrix(x)
  if (is.null(rownames(values)) || any(!nzchar(rownames(values))) ||
      anyDuplicated(rownames(values))) {
    stop("An in-memory character matrix needs unique, non-empty row names.", call. = FALSE)
  }
  blank <- !is.na(values) & !nzchar(trimws(as.character(values)))
  numeric_values <- suppressWarnings(matrix(
    as.numeric(values),
    nrow = nrow(values),
    ncol = ncol(values),
    dimnames = dimnames(values)
  ))
  invalid <- is.na(numeric_values) & !is.na(values) & !blank
  if (any(invalid) || any(numeric_values != round(numeric_values), na.rm = TRUE)) {
    stop("Character states must be integers, NA, or blank.", call. = FALSE)
  }
  numeric_values[is.na(numeric_values)] <- -1
  storage.mode(numeric_values) <- "integer"
  numeric_values
}

.interactive_reference_tree <- function(x) {
  if (inherits(x, "phylo")) return(x)
  if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
    stop("reference_tree must be a Newick path or an ape::phylo object.", call. = FALSE)
  }
  tree <- ape::read.tree(x)
  if (inherits(tree, "multiPhylo") || !inherits(tree, "phylo")) {
    stop("reference_tree must contain exactly one tree.", call. = FALSE)
  }
  tree
}

.interactive_method_args <- function(args) {
  defaults <- list(
    iqtree = list(seed = 1L),
    vine = list(nj_only = TRUE),
    cassiopeia = list(solver = "greedy")
  )
  if (is.null(args)) return(defaults)
  if (!is.list(args) || (length(args) && is.null(names(args)))) {
    stop("args must be a named list keyed by tree-building method.", call. = FALSE)
  }
  for (method in names(args)) {
    if (method %in% names(defaults) && is.list(args[[method]])) {
      defaults[[method]] <- utils::modifyList(defaults[[method]], args[[method]])
    } else {
      defaults[[method]] <- args[[method]]
    }
  }
  defaults
}

#' Plot trees returned by run_tree_building_example()
#'
#' @param result Result from run_tree_building_example().
#' @param cex Tip-label size.
#' @return result, invisibly.
plot_tree_building_result <- function(result, cex = 0.8) {
  if (!is.list(result) || !inherits(result$reference, "phylo") ||
      !is.list(result$trees)) {
    stop("result must come from run_tree_building_example().", call. = FALSE)
  }
  trees <- result$trees[!vapply(result$trees, is.null, logical(1L))]
  candidates <- c(list(reference = result$reference), trees)
  candidates <- lapply(candidates, .tree_building_helpers$first_tree)
  columns <- min(2L, length(candidates))
  rows <- ceiling(length(candidates) / columns)
  old_par <- graphics::par(c("mfrow", "mar"))
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(rows, columns), mar = c(1, 1, 3, 1))
  for (method in names(candidates)) {
    ape::plot.phylo(candidates[[method]], cex = cex, no.margin = TRUE)
    graphics::title(main = method)
  }
  invisible(result)
}

#' Build and compare lineage trees interactively
#'
#' @param reference_tree Path to one Newick tree, or an ape::phylo object.
#' @param character_matrix Path to a TSV character matrix, or a matrix/data
#'   frame with cells as rows and targets as columns.
#' @param methods cliqueR reconstruction methods passed to build_trees().
#' @param output_dir Optional directory for Newick, TSV, and PDF output. NULL
#'   keeps the analysis in memory.
#' @param args Named per-method arguments passed to build_trees().
#' @param on_error Whether a failed optional backend warns or stops.
#' @param plot Whether to draw the reference and reconstructed trees.
#' @return Invisibly, a list containing reference, characters, trees,
#'   comparison, pairwise_rf, and skipped_methods.
run_tree_building_example <- function(
    reference_tree = tree_building_example_files()[["reference_tree"]],
    character_matrix = tree_building_example_files()[["character_matrix"]],
    methods = c("nj", "parsimony"),
    output_dir = NULL,
    args = NULL,
    on_error = c("warn", "stop"),
    plot = interactive()) {
  if (!requireNamespace("cliqueR", quietly = TRUE)) {
    stop("Install cliqueR before running this example.", call. = FALSE)
  }
  for (package in c("ape", "phangorn")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the suggested R package '", package, "'.", call. = FALSE)
    }
  }
  on_error <- match.arg(on_error)

  reference <- .interactive_reference_tree(reference_tree)
  characters <- .interactive_character_matrix(character_matrix)
  if (anyDuplicated(reference$tip.label)) {
    stop("Reference-tree tip labels must be unique.", call. = FALSE)
  }
  missing_tips <- setdiff(reference$tip.label, rownames(characters))
  extra_cells <- setdiff(rownames(characters), reference$tip.label)
  if (length(missing_tips) || length(extra_cells)) {
    stop(
      "Reference tips and matrix cells differ. Missing from matrix: ",
      paste(missing_tips, collapse = ", "),
      "; absent from tree: ", paste(extra_cells, collapse = ", "),
      call. = FALSE
    )
  }
  characters <- characters[reference$tip.label, , drop = FALSE]

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ape::write.tree(reference, file.path(output_dir, "tree_reference.nwk"))
  }
  trees <- cliqueR::build_trees(
    characters,
    methods = methods,
    output_dir = output_dir,
    args = .interactive_method_args(args),
    on_error = on_error
  )
  successful <- trees[!vapply(trees, is.null, logical(1L))]
  if (length(successful) == 0L) {
    stop("No tree-building method succeeded.", call. = FALSE)
  }

  comparison <- .tree_building_helpers$compare_to_reference(
    reference,
    successful,
    characters
  )
  pairwise <- .tree_building_helpers$pairwise_rf(reference, successful)
  result <- list(
    reference = reference,
    characters = characters,
    trees = trees,
    comparison = comparison,
    pairwise_rf = pairwise,
    skipped_methods = setdiff(methods, names(successful))
  )
  class(result) <- c("cliqueR_tree_example", "list")

  if (!is.null(output_dir)) {
    utils::write.table(
      comparison,
      file.path(output_dir, "tree_comparison.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(
      data.frame(tree = rownames(pairwise), pairwise, check.names = FALSE),
      file.path(output_dir, "pairwise_normalized_rf.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(
      data.frame(cell_id = rownames(characters), characters, check.names = FALSE),
      file.path(output_dir, "character_matrix.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    .tree_building_helpers$plot_trees(
      reference,
      successful,
      file.path(output_dir, "tree_comparison.pdf")
    )
  }
  if (isTRUE(plot)) plot_tree_building_result(result)

  cat("\nTree comparison (lower RF and parsimony; higher correlation are better):\n")
  print(comparison, row.names = FALSE, digits = 3)
  if (length(result$skipped_methods)) {
    cat(
      "\nSkipped methods: ",
      paste(result$skipped_methods, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  invisible(result)
}

message(
  "Loaded run_tree_building_example(); run ",
  "result <- run_tree_building_example() to start."
)
