test_that("collapse() has the expected signature", {
  expect_args(collapse, c(
    "input_bam", "read_structure", "output_bam",
    "temp_dir", "max_deletion", "find_inversions",
    "fast_reference_lookup", "correct_only", "threads"
  ))
})

test_that("collapse() runs end-to-end on a fixture BAM", {
  skip("implementation pending — needs a small bundled aligned BAM")
})

test_that("collapse() boolean flags only appear when TRUE", {
  skip("implementation pending — find_inversions=FALSE should not produce --find_inversions")
})

test_that("collapse() correct_only short-circuits consensus building", {
  skip("implementation pending — verify output BAM exists but contains no consensus reads")
})

test_that("collapse() respects a user-supplied temp_dir", {
  skip("implementation pending")
})
