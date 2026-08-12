test_that("tree_nj() has the expected signature", {
  expect_args(tree_nj, c("d", "method"))
})

make_dist <- function() {
  m <- rbind(
    a = c(0L, 1L, 2L),
    b = c(0L, 1L, 0L),
    c = c(1L, 0L, 0L),
    d = c(2L, 2L, 2L)
  )
  indel_distance(m)
}

test_that("tree_nj() returns a phylo with all input cells as tips", {
  skip_if_not_installed("ape")
  tr <- tree_nj(make_dist())
  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, c("a", "b", "c", "d"))
})

test_that("tree_nj() preserves cell barcodes as tip labels", {
  skip_if_not_installed("ape")
  d <- make_dist()
  tr <- tree_nj(d)
  expect_setequal(tr$tip.label, labels(d))
})

test_that("tree_nj() dispatches to bionj when method='bionj'", {
  skip_if_not_installed("ape")
  d <- make_dist()
  expect_identical(tree_nj(d, "bionj")$edge, ape::bionj(d)$edge)
})

test_that("tree_nj() errors when given a non-dist input", {
  expect_error(tree_nj(matrix(0, 3, 3)), "dist")
})
