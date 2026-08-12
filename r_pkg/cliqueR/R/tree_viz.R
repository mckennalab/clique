# Validated colorblind-safe categorical palette (dataviz skill, light surface):
# passes the CVD-separation and lightness-band checks; light hues are relieved
# with a thin tile border and a legend so identity is never colour-alone.
.clique_categorical_palette <- c(
  "#2a78d6", "#1baf7a", "#eda100", "#008300",
  "#4a3aa7", "#e34948", "#e87ba4", "#eb6834"
)
.clique_surface <- "#fcfcfb"   # chart surface (tile borders)
.clique_branch  <- "#52514e"   # secondary ink (branches)
.clique_ink      <- "#0b0b0b"  # primary ink
.clique_ink2     <- "#52514e"  # secondary ink (text)
.clique_na_fill  <- "grey85"

#' The cliqueR default categorical colour palette
#'
#' Eight colourblind-safe hues (validated for colour-vision-deficiency separation
#' on a light surface), assigned to categories in order.
#' @return A character vector of 8 hex colours.
#' @export
clique_categorical_palette <- function() .clique_categorical_palette

#' @noRd
.require_viz_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    rlang::abort(c(
      sprintf("plot_lineage_tree() needs package(s) that are not installed: %s",
              paste(missing, collapse = ", ")),
      i = "ggtree / ggtreeExtra are Bioconductor packages: BiocManager::install(c('ggtree','ggtreeExtra'))",
      i = "ggnewscale is on CRAN: install.packages('ggnewscale')"
    ))
  }
}

#' Coerce `tree` (Newick path, Newick string, phylo, or treedata) to a tree.
#' @noRd
.as_phylo <- function(tree) {
  if (inherits(tree, c("phylo", "treedata"))) return(tree)
  if (is.character(tree) && length(tree) == 1L) {
    if (file.exists(tree)) return(ape::read.tree(tree))
    return(ape::read.tree(text = tree))
  }
  rlang::abort("`tree` must be a `phylo`/`treedata` object, a Newick file path, or a Newick string.")
}

#' Pick the annotation column whose values best match the tree tip labels.
#' @noRd
.auto_mapping_column <- function(ann, tips) {
  overlaps <- vapply(names(ann), function(cn) sum(as.character(ann[[cn]]) %in% tips), integer(1))
  if (length(overlaps) == 0 || max(overlaps) == 0) {
    rlang::abort(c(
      "Could not auto-detect a tip-mapping column: no annotation column matches the tree tip labels.",
      i = "Pass `mapping_column` naming the column that holds the tip identifiers."
    ))
  }
  best <- names(ann)[which.max(overlaps)]
  rlang::inform(sprintf("Mapping annotations to tips by column '%s' (%d/%d tips matched).",
                        best, max(overlaps), length(tips)))
  best
}

#' Named level -> colour map for a categorical column.
#' @noRd
.categorical_colors <- function(levels, palette) {
  n <- length(levels)
  if (n <= length(palette)) {
    cols <- palette[seq_len(n)]
  } else {
    rlang::warn(sprintf(
      "%d categories exceed the %d colourblind-safe palette colours; extending with an HCL qualitative palette (less distinguishable). Consider grouping rare levels.",
      n, length(palette)))
    cols <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  stats::setNames(cols, levels)
}

#' Continuous fill scale from a viridis option name or a vector of gradient colours.
#' @noRd
.continuous_scale <- function(name, cont) {
  if (length(cont) == 1L) {
    ggplot2::scale_fill_viridis_c(name = name, option = cont, na.value = .clique_na_fill)
  } else {
    ggplot2::scale_fill_gradientn(name = name, colours = cont, na.value = .clique_na_fill)
  }
}

#' Plot a lineage tree with per-tip annotation boxes
#'
#' Renders a (by default circular) lineage tree with a concentric ring of filled
#' boxes at the leaves for each annotation column. Categorical annotations are
#' coloured with a colourblind-safe discrete palette; continuous annotations use
#' a gradient. Each annotation gets its own fill scale and legend.
#'
#' @param tree A tree: an `ape::phylo` / `treedata` object, a path to a Newick
#'   file, or a Newick string.
#' @param annotations A data frame of per-tip annotations (one row per tip),
#'   or `NULL` to draw the bare tree.
#' @param columns Which annotation columns to draw as rings, innermost first.
#'   Defaults to every column except `mapping_column`.
#' @param mapping_column The `annotations` column holding the tip identifiers
#'   (matched against the tree's tip labels). Auto-detected when `NULL`.
#' @param layout One of `"circular"` (default), `"fan"`, or `"rectangular"`.
#' @param categorical_palette Colours for categorical annotations, assigned to
#'   levels in order. Defaults to [clique_categorical_palette()].
#' @param continuous_palette Gradient for continuous annotations: a viridis
#'   option name (e.g. `"viridis"`, `"magma"`) or a vector of >= 2 colours.
#' @param ring_width Radial width of each annotation ring (fraction of the tree).
#' @param ring_gap Gap before each ring (fraction of the tree).
#' @param tip_labels Draw tip labels?
#' @param title Optional plot title.
#' @return A `ggplot`/`ggtree` object (print it, or save with [ggplot2::ggsave()]).
#' @examples
#' \dontrun{
#' plot_lineage_tree(
#'   "tree.nwk",
#'   annotations    = cell_meta,       # data frame with a tip-id column
#'   mapping_column = "cell_barcode",
#'   columns        = c("cluster", "expression"),  # categorical + continuous
#'   layout         = "circular"
#' )
#' }
#' @export
plot_lineage_tree <- function(tree,
                              annotations = NULL,
                              columns = NULL,
                              mapping_column = NULL,
                              layout = c("circular", "fan", "rectangular"),
                              categorical_palette = clique_categorical_palette(),
                              continuous_palette = "viridis",
                              ring_width = 0.12,
                              ring_gap = 0.04,
                              tip_labels = FALSE,
                              title = NULL) {
  .require_viz_pkgs(c("ggtree", "ggplot2"))
  layout <- match.arg(layout)
  phy <- .as_phylo(tree)
  tips <- if (inherits(phy, "treedata")) phy@phylo$tip.label else phy$tip.label

  p <- ggtree::ggtree(phy, layout = layout, linewidth = 0.35, color = .clique_branch)
  if (isTRUE(tip_labels)) {
    p <- p + ggtree::geom_tiplab(size = 2, color = .clique_ink2, offset = ring_gap)
  }

  if (!is.null(annotations)) {
    .require_viz_pkgs("ggnewscale")
    annotations <- as.data.frame(annotations, stringsAsFactors = FALSE)

    if (is.null(mapping_column)) mapping_column <- .auto_mapping_column(annotations, tips)
    if (!mapping_column %in% names(annotations)) {
      rlang::abort(sprintf("`mapping_column` '%s' is not a column of `annotations`.", mapping_column))
    }
    ids <- as.character(annotations[[mapping_column]])
    if (anyDuplicated(ids)) {
      rlang::abort(sprintf("`mapping_column` '%s' has duplicate values; tip identifiers must be unique.", mapping_column))
    }

    if (length(setdiff(ids, tips))) {
      rlang::warn(sprintf("%d annotation row(s) match no tree tip and are ignored (e.g. %s).",
                          length(setdiff(ids, tips)),
                          paste(utils::head(setdiff(ids, tips), 3), collapse = ", ")))
    }
    if (length(setdiff(tips, ids))) {
      rlang::warn(sprintf("%d tip(s) have no annotation and are drawn in the NA colour.",
                          length(setdiff(tips, ids))))
    }

    if (is.null(columns)) columns <- setdiff(names(annotations), mapping_column)
    bad <- setdiff(columns, names(annotations))
    if (length(bad)) rlang::abort(sprintf("Column(s) not in `annotations`: %s", paste(bad, collapse = ", ")))
    if (length(columns) == 0) rlang::abort("No annotation columns to draw.")

    # Each annotation is one concentric ring of tiles (via gheatmap), given its
    # own fill scale + legend with ggnewscale. Offsets are in the tree's x units.
    xspan <- diff(range(p$data$x, na.rm = TRUE))
    if (!is.finite(xspan) || xspan == 0) xspan <- 1
    ring_x <- ring_width * xspan
    gap_x  <- ring_gap * xspan

    for (i in seq_along(columns)) {
      col <- columns[i]
      vals <- annotations[[col]]
      is_cat <- !(is.numeric(vals) && !is.factor(vals))
      d1 <- data.frame(value = if (is_cat) factor(vals) else vals, row.names = ids)
      names(d1) <- col
      d1 <- d1[rownames(d1) %in% tips, , drop = FALSE]

      # Build the fill scale first, so a "too many categories" warning is not
      # swallowed by the message/warning suppression around gheatmap's internals.
      fill_scale <- if (is_cat) {
        ggplot2::scale_fill_manual(
          name = col,
          values = .categorical_colors(levels(d1[[col]]), categorical_palette),
          na.value = .clique_na_fill, drop = FALSE
        )
      } else {
        .continuous_scale(col, continuous_palette)
      }

      # gheatmap emits benign ggtree-internal messages/warnings; mute those here.
      if (i > 1) p <- suppressMessages(p + ggnewscale::new_scale_fill())
      p <- suppressMessages(suppressWarnings(ggtree::gheatmap(
        p, d1,
        offset   = gap_x + (i - 1) * (ring_x + gap_x),
        width    = ring_width,
        color    = .clique_surface,   # thin tile border -> relief for light hues
        colnames = FALSE              # the legend names each ring
      )))
      p <- suppressMessages(p + fill_scale)
    }
  }

  p <- p + ggplot2::theme(
    legend.title = ggplot2::element_text(size = 9, color = .clique_ink),
    legend.text  = ggplot2::element_text(size = 8, color = .clique_ink2),
    plot.title   = ggplot2::element_text(size = 12, color = .clique_ink)
  )
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}
