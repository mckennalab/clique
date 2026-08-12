test_that("align() has the expected signature", {
  expect_args(align, c(
    "read1", "read2", "index1", "index2",
    "read_structure", "output_bam",
    "aligner", "max_reference_multiplier",
    "min_read_length", "threads"
  ))
})

test_that("align() runs end-to-end on a fixture FASTQ", {
  skip("implementation pending — needs a small bundled FASTQ + read-structure YAML")
})

test_that("align() drops NULL paired/index reads from the argv", {
  skip("implementation pending — verify single-end call does not pass --read2/--index1/--index2")
})

test_that("align() forwards the aligner choice to the CLI", {
  skip("implementation pending — assert --aligner Degenerate appears when aligner='Degenerate'")
})

test_that("align() errors clearly when read1 is missing or unreadable", {
  skip("implementation pending")
})
