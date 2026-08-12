test_that("tree backends validate integer matrices and cell labels", {
  expect_error(cliqueR:::.validate_indel_matrix(matrix(1:4, 2)), "at least 3")

  m <- make_lineage_matrix()
  m[1, 1] <- 0.5
  expect_error(cliqueR:::.validate_indel_matrix(m), "must be integers")

  m <- make_lineage_matrix()
  rownames(m)[2] <- rownames(m)[1]
  expect_error(cliqueR:::.validate_indel_matrix(m), "must be unique")
})

test_that("all-missing sites are removed without changing cell names", {
  m <- cbind(make_lineage_matrix(), empty = -1L)
  expect_warning(valid <- cliqueR:::.validate_indel_matrix(m), "Dropping 1 site")
  expect_false("empty" %in% colnames(valid))
  expect_identical(rownames(valid), rownames(m))
})

test_that("safe external labels round-trip arbitrary cell names", {
  skip_if_not_installed("ape")
  labels <- rownames(make_lineage_matrix())
  tip_map <- cliqueR:::.tree_tip_map(labels)
  tr <- ape::read.tree(text = safe_four_tip_newick())
  restored <- cliqueR:::.restore_tree_tips(tr, tip_map, "test")
  expect_setequal(restored$tip.label, labels)
})
