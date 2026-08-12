test_that("tree_mix() has the expected signature", {
  expect_args(
    tree_mix,
    c("indel_mat", "method", "outgroup", "all", "binary",
      "output_prefix", "missing", "timeout")
  )
})

test_that("tree_mix() one-hot encodes states and restores labels", {
  skip_on_os("windows")
  fake <- make_fake_executable(c(
    "cat > menu.captured",
    sprintf("printf '%%s\\n%%s\\n' '%s' '%s' > outtree", safe_four_tip_newick(), safe_four_tip_newick()),
    "printf 'fake MIX report\\n' > outfile"
  ))
  prefix <- tempfile("mix-result-")
  tr <- tree_mix(
    make_lineage_matrix(), binary = fake, output_prefix = prefix,
    outgroup = "cell/C"
  )

  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, rownames(make_lineage_matrix()))
  expect_identical(attr(tr, "clique_tied_trees"), 2L)
  input <- readLines(paste0(prefix, ".mix.infile"))
  expect_match(input[1], "4\\s+4")
  expect_true(any(grepl("?", input, fixed = TRUE)))
  menu <- readLines(paste0(prefix, ".mix.menu"))
  expect_identical(menu, c("P", "O", "3", "Y"))
})

test_that("tree_mix(all=TRUE) returns every tied tree", {
  skip_on_os("windows")
  fake <- make_fake_executable(c(
    "cat >/dev/null",
    sprintf("printf '%%s\\n%%s\\n' '%s' '%s' > outtree", safe_four_tip_newick(), safe_four_tip_newick())
  ))
  trees <- tree_mix(make_lineage_matrix(), binary = fake, method = "wagner", all = TRUE)
  expect_s3_class(trees, "multiPhylo")
  expect_length(trees, 2L)
  expect_true(all(vapply(trees, function(x) setequal(x$tip.label, rownames(make_lineage_matrix())), logical(1))))
})

test_that("tree_mix() rejects a matrix without edits", {
  m <- matrix(0L, 4, 3, dimnames = list(letters[1:4], paste0("s", 1:3)))
  fake <- make_fake_executable("exit 0")
  expect_error(tree_mix(m, binary = fake), "positive edit state")
})
