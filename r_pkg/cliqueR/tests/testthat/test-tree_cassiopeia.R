test_that("tree_cassiopeia() has the expected signature", {
  expect_args(tree_cassiopeia, c("indel_mat", "solver", "priors", "missing", "..."))
})

test_that("tree_cassiopeia() errors helpfully when reticulate or cassiopeia is missing", {
  skip_if_not_installed("reticulate")
  skip_if(cassiopeia_test_available(), "Cassiopeia is installed in the active Python environment")
  expect_error(
    suppressWarnings(tree_cassiopeia(make_lineage_matrix())),
    "Python package 'cassiopeia'"
  )
})

test_that("tree_cassiopeia() greedy solver returns a phylo with all cells as tips", {
  skip_if_not_installed("ape")
  skip_if(!cassiopeia_test_available(), "Cassiopeia is not installed")
  tr <- tree_cassiopeia(make_lineage_matrix(), solver = "greedy")
  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, rownames(make_lineage_matrix()))
})

test_that("tree_cassiopeia() passes priors through to the solver", {
  skip_if_not_installed("ape")
  skip_if(!cassiopeia_test_available(), "Cassiopeia is not installed")
  priors <- list(
    target_1 = c(`1` = 0.5, `2` = 0.5),
    target_2 = c(`1` = 1),
    target_3 = c(`1` = 1)
  )
  expect_s3_class(
    tree_cassiopeia(make_lineage_matrix(), solver = "greedy", priors = priors),
    "phylo"
  )
})

test_that("tree_cassiopeia() validates the solver name before calling Python", {
  expect_error(tree_cassiopeia(make_lineage_matrix(), solver = "not-a-solver"),
               "should be one of")
})
