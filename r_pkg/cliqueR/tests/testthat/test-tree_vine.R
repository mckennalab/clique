test_that("tree_vine() has the expected signature", {
  expect_args(
    tree_vine,
    c("indel_mat", "model_type", "mutation_prior", "dimensionality", "threads",
      "nsamples", "nj_only", "return_posterior", "binary",
      "output_prefix", "extra_args", "missing", "timeout")
  )
})

test_that("tree_vine() returns the mean and optional posterior trees", {
  skip_on_os("windows")
  fake <- make_fake_executable(c(
    "mean=''",
    "shift",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    --mean) mean=$2; shift 2 ;;",
    "    *) shift ;;",
    "  esac",
    "done",
    sprintf("printf '%%s\\n' '%s' > \"$mean\"", safe_four_tip_newick()),
    sprintf("printf '%%s\\n%%s\\n' '%s' '%s'", safe_four_tip_newick(), safe_four_tip_newick())
  ))
  prefix <- tempfile("vine-result-")
  tr <- tree_vine(
    make_lineage_matrix(), binary = fake, output_prefix = prefix,
    nsamples = 2L, return_posterior = TRUE
  )

  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, rownames(make_lineage_matrix()))
  expect_s3_class(attr(tr, "clique_posterior"), "multiPhylo")
  expect_length(attr(tr, "clique_posterior"), 2L)
  matrix_lines <- readLines(paste0(prefix, ".vine.tsv"))
  expect_match(matrix_lines[1], "cell\\ttarget_1")
  expect_true(any(grepl("-1", matrix_lines, fixed = TRUE)))
  expect_true(file.exists(paste0(prefix, ".vine.samples.nwk")))
})

test_that("tree_vine() supports NJ-only output", {
  skip_on_os("windows")
  fake <- make_fake_executable(sprintf("printf '%%s\\n' '%s'", safe_four_tip_newick()))
  tr <- tree_vine(make_lineage_matrix(), binary = fake, nj_only = TRUE)
  expect_s3_class(tr, "phylo")
  expect_null(attr(tr, "clique_posterior"))
})

test_that("tree_vine() enforces the native missing value", {
  expect_error(tree_vine(make_lineage_matrix(), missing = -9L), "requires `missing = -1`")
})

test_that("tree_vine() validates the embedding dimension", {
  expect_error(
    tree_vine(make_lineage_matrix(), dimensionality = 4L),
    "smaller than the number of cells"
  )
})
