test_that("read_10x_lineage() has the expected signature", {
  expect_args(read_10x_lineage, c("cellranger_dir", "lineage_bam", "min_reads_per_cell"))
})

test_that("read_10x_lineage() returns a tibble with the documented columns", {
  skip("implementation pending — expect cell_barcode, umi, consensus_seq, n_reads")
})

test_that("read_10x_lineage() drops cells below min_reads_per_cell", {
  skip("implementation pending")
})

test_that("read_10x_lineage() restricts to barcodes in filtered_feature_bc_matrix", {
  skip("implementation pending — reads tagged with barcodes not in the matrix should be ignored")
})

test_that("read_10x_lineage() errors clearly when cellranger_dir is malformed", {
  skip("implementation pending")
})
