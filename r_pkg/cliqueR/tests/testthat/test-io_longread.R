test_that("read_longread_lineage() has the expected signature", {
  expect_args(read_longread_lineage,
              c("bam_path", "read_structure", "tags", "region", "samtools"))
})

# Build a tiny two-read BAM fixture from a SAM string, using samtools.
make_fixture_bam <- function() {
  skip_if(unname(Sys.which("samtools")) == "", "samtools not on PATH")
  sam <- c(
    "@HD\tVN:1.6\tSO:unsorted",
    "@SQ\tSN:ref\tLN:20",
    paste("r1", 0, "ref", 1, 255, "20M", "*", 0, 0,
          "ACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIII",
          "ce:Z:NONE", "rc:i:1", sep = "\t"),
    paste("r2", 0, "ref", 1, 255, "20M", "*", 0, 0,
          "ACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIII",
          "ce:Z:2D+10", "rc:i:3", sep = "\t")
  )
  sam_path <- tempfile(fileext = ".sam")
  bam_path <- tempfile(fileext = ".bam")
  writeLines(sam, sam_path)
  processx::run("samtools", c("view", "-b", "-o", bam_path, sam_path))
  bam_path
}

test_that("read_longread_lineage() returns a data frame keyed by read name", {
  bam <- make_fixture_bam()
  df <- read_longread_lineage(bam, tags = c("ce", "rc"))
  expect_identical(nrow(df), 2L)
  expect_false(any(duplicated(df$read_name)))
  expect_identical(df$read_name, c("r1", "r2"))
  expect_identical(unique(df$reference), "ref")
})

test_that("read_longread_lineage() extracts requested BAM tags as columns", {
  bam <- make_fixture_bam()
  df <- read_longread_lineage(bam, tags = c("ce", "rc"))
  expect_true(all(c("ce", "rc") %in% names(df)))
  expect_identical(df$ce, c("NONE", "2D+10"))
  expect_identical(df$rc, c("1", "3"))
})

test_that("read_longread_lineage() returns all parsed tags when tags=NULL", {
  bam <- make_fixture_bam()
  df <- read_longread_lineage(bam, tags = NULL)
  expect_true(all(c("ce", "rc") %in% names(df)))
})

test_that("read_longread_lineage() errors on a missing BAM file", {
  expect_error(read_longread_lineage("/no/such/file.bam"), "not found")
})
