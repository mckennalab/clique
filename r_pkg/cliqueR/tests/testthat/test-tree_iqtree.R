test_that("tree_iqtree() has the expected signature", {
  expect_args(
    tree_iqtree,
    c("indel_mat", "model", "threads", "seed", "outgroup", "binary",
      "output_prefix", "extra_args", "missing", "timeout")
  )
})

test_that("tree_iqtree() writes MORPH NEXUS and restores labels", {
  skip_on_os("windows")
  fake <- make_fake_executable(c(
    "prefix=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    -pre) prefix=$2; shift 2 ;;",
    "    *) shift ;;",
    "  esac",
    "done",
    sprintf("printf '%%s\\n' '%s' > \"${prefix}.treefile\"", safe_four_tip_newick())
  ))
  prefix <- tempfile("iqtree-result-")
  tr <- tree_iqtree(
    make_lineage_matrix(), binary = fake, output_prefix = prefix,
    outgroup = "cell:A"
  )

  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, rownames(make_lineage_matrix()))
  expect_identical(attr(tr, "clique_backend"), "iqtree2")
  nex <- readLines(paste0(prefix, ".nex"))
  expect_true(any(grepl("DATATYPE=STANDARD", nex, fixed = TRUE)))
  expect_true(any(grepl("MISSING=?", nex, fixed = TRUE)))
  expect_true(any(grepl("clq0000001", nex, fixed = TRUE)))
  expect_true(any(grepl("?", nex, fixed = TRUE)))
})

test_that("tree_iqtree() rejects sites with more than 32 states", {
  m <- matrix(0:32, ncol = 1, dimnames = list(paste0("c", 0:32), "target"))
  fake <- make_fake_executable("exit 0")
  expect_error(tree_iqtree(m, binary = fake), "at most 32")
})

test_that("tree_iqtree() protects wrapper-managed arguments", {
  expect_error(
    tree_iqtree(make_lineage_matrix(), binary = "/missing", extra_args = c("-pre", "x")),
    "cannot override"
  )
})
