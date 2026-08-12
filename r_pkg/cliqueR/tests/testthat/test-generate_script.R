test_that("generate_clique_script() has the expected signature", {
  expect_args(generate_clique_script,
              c("read_structure", "reference", "read1", "output_script",
                "read2", "index1", "index2", "output_dir", "sample_name",
                "clique_bin", "samtools", "threads", "aligner",
                "min_read_length", "max_reference_multiplier", "collapse",
                "correct_only", "sort_and_index", "validate_reference", "execute"))
})

# Minimal fixtures: a read-structure YAML, a matching FASTA, and dummy read files.
make_fixtures <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  yaml <- file.path(d, "rs.yaml")
  writeLines(c(
    "---", "known_strand: true", "reads:", "  - !Read1", "    orientation: Forward",
    "references:",
    "  refA:", "    sequence: \"ACGTACGT\"", "    targets: []", "    target_types: []",
    "    umi_configurations: {}",
    "  refB:", "    sequence: \"TTTTGGGG\"", "    targets: []", "    target_types: []",
    "    umi_configurations: {}"
  ), yaml)
  fa <- file.path(d, "ref.fa")
  writeLines(c(">refA", "ACGTACGT", ">refB", "TTTTGGGG"), fa)
  r1 <- file.path(d, "s_R1.fastq.gz"); file.create(r1)
  r2 <- file.path(d, "s_R2.fastq.gz"); file.create(r2)
  list(dir = d, yaml = yaml, fa = fa, r1 = r1, r2 = r2)
}

test_that("the YAML reference-name parser finds the reference keys", {
  f <- make_fixtures()
  expect_setequal(cliqueR:::.yaml_reference_names(f$yaml), c("refA", "refB"))
})

test_that("a basic script contains align, sort/index and the correct flags", {
  f <- make_fixtures()
  out <- file.path(f$dir, "run.sh")
  suppressMessages(generate_clique_script(
    read_structure = f$yaml, reference = f$fa, read1 = f$r1,
    output_script = out, output_dir = f$dir, clique_bin = "clique", threads = 3
  ))
  s <- paste(readLines(out), collapse = "\n")
  expect_match(s, 'CLIQUE" align')
  expect_match(s, "--read-structure")
  expect_match(s, "--aligner wfa")          # lower-cased for clap
  expect_match(s, "--threads 3")
  expect_match(s, 'SAMTOOLS" sort')
  expect_match(s, 'SAMTOOLS" index')
  expect_no_match(s, 'CLIQUE" collapse')    # collapse defaults to FALSE
  # the script is marked executable
  expect_true(file.access(out, mode = 1L) == 0L)
})

test_that("optional read/index flags appear only when supplied", {
  f <- make_fixtures()
  out <- file.path(f$dir, "run.sh")
  suppressMessages(generate_clique_script(
    read_structure = f$yaml, reference = f$fa, read1 = f$r1, read2 = f$r2,
    output_script = out, output_dir = f$dir, clique_bin = "clique"
  ))
  s <- paste(readLines(out), collapse = "\n")
  expect_match(s, "--read2")
  expect_no_match(s, "--index1")
})

test_that("collapse adds the collapse step and forces sort/index", {
  f <- make_fixtures()
  out <- file.path(f$dir, "run.sh")
  suppressMessages(generate_clique_script(
    read_structure = f$yaml, reference = f$fa, read1 = f$r1, output_script = out,
    output_dir = f$dir, clique_bin = "clique", collapse = TRUE, correct_only = TRUE,
    sort_and_index = FALSE
  ))
  s <- paste(readLines(out), collapse = "\n")
  expect_match(s, 'CLIQUE" collapse')
  expect_match(s, "--correct-only")
  expect_match(s, 'SAMTOOLS" index')         # forced on despite sort_and_index = FALSE
})

test_that("mismatched reference names produce a warning", {
  f <- make_fixtures()
  out <- file.path(f$dir, "run.sh")
  bad_fa <- file.path(f$dir, "bad.fa")
  writeLines(c(">refA", "ACGTACGT", ">refZ", "GGGGCCCC"), bad_fa)  # refZ not in YAML
  expect_warning(
    suppressMessages(generate_clique_script(
      read_structure = f$yaml, reference = bad_fa, read1 = f$r1,
      output_script = out, output_dir = f$dir, clique_bin = "clique"
    )),
    "Reference names differ"
  )
})

test_that("missing inputs and bad aligner error clearly", {
  f <- make_fixtures()
  out <- file.path(f$dir, "run.sh")
  expect_error(
    generate_clique_script(f$yaml, f$fa, file.path(f$dir, "nope.fastq.gz"), out,
                           clique_bin = "clique"),
    "read1 FASTQ not found"
  )
  expect_error(
    suppressMessages(generate_clique_script(f$yaml, f$fa, f$r1, out,
                     clique_bin = "clique", aligner = "bowtie")),
    "Unknown aligner"
  )
})
