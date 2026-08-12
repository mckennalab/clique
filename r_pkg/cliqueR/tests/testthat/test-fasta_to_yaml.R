test_that("fasta_targets_to_yaml() has the expected signature", {
  expect_args(fasta_targets_to_yaml,
              c("fasta_file", "targets", "output_file", "target_types",
                "search_reverse_complement", "keep_references_without_targets",
                "merge", "known_strand", "reads", "uppercase"))
})

# A two-record FASTA: ref_a contains target T1 forward; ref_b contains the
# reverse complement of T2.
write_fixture_fasta <- function() {
  path <- tempfile(fileext = ".fa")
  writeLines(c(
    ">ref_a some description",
    "AAAAAGACGGCTATACAAGGCATCGCGGTTTTT",
    ">ref_b",
    "CCCCCAAAAAAAAAAAAAAAAAAAAAAAAAGGGGG"
  ), path)
  path
}

test_that("forward targets are placed in the references that contain them", {
  fa <- write_fixture_fasta()
  out <- tempfile(fileext = ".yaml")
  res <- suppressMessages(fasta_targets_to_yaml(
    fasta_file = fa,
    targets = "GACGGCTATACAAGGCATCGCGG",
    output_file = out,
    target_types = "Cas9WT"
  ))
  # ref_a matched forward; ref_b kept but empty
  expect_equal(res$ref_a$target, "GACGGCTATACAAGGCATCGCGG")
  expect_equal(res$ref_a$strand, "+")
  expect_length(res$ref_b$target, 0)

  yaml <- readLines(out)
  expect_true(any(grepl("ref_a:", yaml)))
  expect_true(any(grepl('targets: \\["GACGGCTATACAAGGCATCGCGG"\\]', yaml)))
  expect_true(any(grepl("target_types: \\[\\]", yaml)))  # ref_b empty
})

test_that("a reverse-complement-only target is written as its forward substring", {
  fa <- write_fixture_fasta()
  out <- tempfile(fileext = ".yaml")
  # supply the RC of a run that occurs forward in ref_b, so it is only found via RC
  rc_target <- cliqueR:::.revcomp("AAAAAAAAAAAAAAAAAAAAAAAAA")  # poly-T RC of poly-A
  res <- suppressMessages(fasta_targets_to_yaml(
    fasta_file = fa,
    targets = rc_target,
    output_file = out,
    target_types = "Cas9WT"
  ))
  # found in ref_b via reverse complement; the written target is the forward substring
  expect_equal(res$ref_b$strand, "-")
  expect_true(grepl(res$ref_b$target, "CCCCCAAAAAAAAAAAAAAAAAAAAAAAAAGGGGG", fixed = TRUE))
})

test_that("keep_references_without_targets = FALSE drops unmatched references", {
  fa <- write_fixture_fasta()
  out <- tempfile(fileext = ".yaml")
  res <- suppressMessages(fasta_targets_to_yaml(
    fasta_file = fa,
    targets = "GACGGCTATACAAGGCATCGCGG",
    output_file = out,
    keep_references_without_targets = FALSE
  ))
  expect_named(res, "ref_a")
  expect_false("ref_b" %in% names(res))
})

test_that("invalid target types are rejected", {
  fa <- write_fixture_fasta()
  out <- tempfile(fileext = ".yaml")
  expect_error(
    fasta_targets_to_yaml(fa, "GACGGCTATACAAGGCATCGCGG", out, target_types = "NotACas"),
    "Unknown clique target type"
  )
})

test_that("target_types must be length 1 or match targets", {
  fa <- write_fixture_fasta()
  out <- tempfile(fileext = ".yaml")
  expect_error(
    fasta_targets_to_yaml(fa, c("AAAA", "CCCC"), out, target_types = c("Cas9WT", "Cas9ABE", "Cas9CBE")),
    "must be 1 or equal"
  )
})

test_that("a data.frame of targets with a type column is honored", {
  fa <- write_fixture_fasta()
  out <- tempfile(fileext = ".yaml")
  tdf <- data.frame(target = "GACGGCTATACAAGGCATCGCGG", type = "Cas9ABE", stringsAsFactors = FALSE)
  res <- suppressMessages(fasta_targets_to_yaml(fa, tdf, out))
  expect_equal(res$ref_a$type, "Cas9ABE")
})
