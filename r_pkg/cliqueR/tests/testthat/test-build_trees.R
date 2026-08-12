test_that("build_trees() has the expected signature", {
  expect_args(build_trees,
              c("indel_mat", "methods", "output_dir", "args", "missing", "on_error"))
})

test_that("build_trees() runs the dependency-free methods and names the result", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  trees <- build_trees(make_lineage_matrix(), methods = c("nj", "parsimony"))
  expect_named(trees, c("nj", "parsimony"))
  expect_s3_class(trees$nj, "phylo")
  expect_s3_class(trees$parsimony, "phylo")
  expect_setequal(trees$nj$tip.label, rownames(make_lineage_matrix()))
})

test_that("build_trees() passes per-method args (nj -> bionj)", {
  skip_if_not_installed("ape")
  trees <- build_trees(make_lineage_matrix(), methods = "nj",
                       args = list(nj = list(method = "bionj")))
  expect_s3_class(trees$nj, "phylo")
})

test_that("build_trees() writes one Newick per successful method", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")
  dir <- withr::local_tempdir()
  build_trees(make_lineage_matrix(), methods = c("nj", "parsimony"), output_dir = dir)
  expect_true(file.exists(file.path(dir, "tree_nj.nwk")))
  expect_true(file.exists(file.path(dir, "tree_parsimony.nwk")))
  # the file is a readable Newick tree
  tr <- ape::read.tree(file.path(dir, "tree_nj.nwk"))
  expect_s3_class(tr, "phylo")
})

test_that("build_trees() skips a failing backend with on_error='warn'", {
  skip_if_not_installed("ape")
  # point iqtree at a non-existent binary so it fails; nj still succeeds
  expect_warning(
    trees <- build_trees(make_lineage_matrix(),
                         methods = c("nj", "iqtree"),
                         args = list(iqtree = list(binary = "/no/such/iqtree"))),
    "iqtree")
  expect_s3_class(trees$nj, "phylo")
  expect_null(trees$iqtree)
})

test_that("build_trees() aborts on a failing method with on_error='stop'", {
  skip_if_not_installed("ape")
  expect_error(
    build_trees(make_lineage_matrix(),
                methods = "iqtree",
                args = list(iqtree = list(binary = "/no/such/iqtree")),
                on_error = "stop"),
    "iqtree")
})

test_that("build_trees() rejects unknown methods and mis-keyed args", {
  expect_error(build_trees(make_lineage_matrix(), methods = "bogus"), "Unknown")
  expect_error(build_trees(make_lineage_matrix(), methods = "nj",
                           args = list(bogus = list())), "method names")
  expect_error(build_trees(make_lineage_matrix(), methods = "nj",
                           args = list(list(method = "bionj"))), "named list")
})
