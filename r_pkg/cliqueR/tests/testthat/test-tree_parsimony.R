test_that("tree_parsimony() has the expected signature", {
  expect_args(tree_parsimony, c("indel_mat", "start", "trace"))
})

test_that("tree_parsimony() returns a phylo with all cells as tips", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  tr <- tree_parsimony(make_lineage_matrix())
  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, rownames(make_lineage_matrix()))
  expect_identical(attr(tr, "clique_backend"), "phangorn-parsimony")
})

test_that("tree_parsimony() accepts a user-supplied start tree", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  m <- make_lineage_matrix()
  start <- tree_nj(indel_distance(m))
  expect_s3_class(tree_parsimony(m, start = start), "phylo")
})

test_that("tree_parsimony() handles missing data (-1) without crashing", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  expect_s3_class(tree_parsimony(make_lineage_matrix()), "phylo")
})

test_that("tree_parsimony() rejects a start tree with different tips", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  bad <- ape::rtree(4)
  expect_error(tree_parsimony(make_lineage_matrix(), start = bad), "same cell names")
})
