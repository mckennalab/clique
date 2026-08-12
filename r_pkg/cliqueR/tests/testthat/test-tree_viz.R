test_that("plot_lineage_tree() has the expected signature", {
  expect_args(plot_lineage_tree,
              c("tree", "annotations", "columns", "mapping_column", "layout",
                "categorical_palette", "continuous_palette", "ring_width",
                "ring_gap", "tip_labels", "title"))
})

test_that("clique_categorical_palette() is 8 valid, unique hex colours", {
  pal <- clique_categorical_palette()
  expect_length(pal, 8)
  expect_true(all(grepl("^#[0-9a-fA-F]{6}$", pal)))
  expect_equal(length(unique(pal)), 8)
})

make_tree <- function() {
  set.seed(1)
  tr <- ape::rtree(8)
  ann <- data.frame(
    id  = tr$tip.label,
    grp = rep(c("x", "y"), length.out = 8),   # categorical
    val = seq(0, 1, length.out = 8),           # continuous
    stringsAsFactors = FALSE
  )
  list(tr = tr, ann = ann)
}

test_that("returns a ggplot object (bare tree and annotated)", {
  skip_if_not_installed("ggtree")
  skip_if_not_installed("ggnewscale")
  f <- make_tree()
  expect_s3_class(plot_lineage_tree(f$tr), "ggplot")
  p <- suppressMessages(plot_lineage_tree(f$tr, annotations = f$ann, mapping_column = "id"))
  expect_s3_class(p, "ggplot")
})

test_that("accepts a Newick string as well as a phylo object", {
  skip_if_not_installed("ggtree")
  f <- make_tree()
  nwk <- ape::write.tree(f$tr)
  expect_s3_class(plot_lineage_tree(nwk), "ggplot")
})

test_that("auto-detects the tip-mapping column", {
  skip_if_not_installed("ggtree")
  skip_if_not_installed("ggnewscale")
  f <- make_tree()
  expect_message(
    plot_lineage_tree(f$tr, annotations = f$ann, columns = "grp"),
    "Mapping annotations to tips by column 'id'"
  )
})

test_that("errors when mapping_column is not a column", {
  skip_if_not_installed("ggtree")
  f <- make_tree()
  expect_error(
    plot_lineage_tree(f$tr, annotations = f$ann, mapping_column = "nope"),
    "not a column"
  )
})

test_that("warns about unannotated tips", {
  skip_if_not_installed("ggtree")
  skip_if_not_installed("ggnewscale")
  f <- make_tree()
  ann2 <- f$ann[1:6, ]   # two tips left unannotated
  expect_warning(
    suppressMessages(plot_lineage_tree(f$tr, annotations = ann2, mapping_column = "id")),
    "no annotation"
  )
})

test_that("warns about annotation rows that are not on the tree", {
  skip_if_not_installed("ggtree")
  skip_if_not_installed("ggnewscale")
  f <- make_tree()
  ann3 <- rbind(f$ann, data.frame(id = "ghost", grp = "z", val = 0.5))
  expect_warning(
    suppressMessages(plot_lineage_tree(f$tr, annotations = ann3, mapping_column = "id")),
    "match no tree tip"
  )
})

test_that("errors on duplicate tip identifiers", {
  skip_if_not_installed("ggtree")
  f <- make_tree()
  ann4 <- f$ann
  ann4$id[2] <- ann4$id[1]   # duplicate id
  expect_error(
    plot_lineage_tree(f$tr, annotations = ann4, mapping_column = "id"),
    "duplicate"
  )
})
