test_that("tree metric functions have the expected signatures", {
  expect_args(rf_distance, c("tree1", "tree2", "normalize"))
  expect_args(triplet_correctness, c("reconstructed", "ground_truth", "n_triplets"))
  expect_args(quartet_distance, c("tree1", "tree2"))
  expect_args(ancestor_recall, c("reconstructed", "ground_truth"))
})

test_that("rf_distance() returns 0 for a tree against itself", {
  skip("implementation pending")
})

test_that("rf_distance() normalize=TRUE bounds output to [0, 1]", {
  skip("implementation pending")
})

test_that("rf_distance() errors when tip-label sets differ", {
  skip("implementation pending — silent intersection would mask data-pipeline bugs")
})

test_that("triplet_correctness() returns 1.0 for a tree against itself", {
  skip("implementation pending")
})

test_that("triplet_correctness() reports n_sampled <= n_triplets", {
  skip("implementation pending")
})

test_that("triplet_correctness() enumerates all triplets when n_triplets=Inf", {
  skip("implementation pending — n_sampled should equal choose(n_tips, 3)")
})

test_that("quartet_distance() returns 0 for a tree against itself", {
  skip("implementation pending")
})

test_that("ancestor_recall() returns 1.0 for a tree against itself", {
  skip("implementation pending")
})

test_that("ancestor_recall() returns 0 when the reconstructed tree is a star", {
  skip("implementation pending — star tree has no informative internal nodes")
})
