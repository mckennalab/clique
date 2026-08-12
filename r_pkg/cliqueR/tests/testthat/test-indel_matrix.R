test_that("build_indel_matrix() has the expected signature", {
  expect_args(build_indel_matrix,
              c("lineage_df", "sites", "cell_col", "event_col",
                "reference_col", "missing_threshold", "missing_value"))
})

test_that("indel_distance() has the expected signature", {
  expect_args(indel_distance, c("indel_mat", "weights", "missing"))
})

# A 5-read, 2-site fixture. site1: NONE/2D+10, site2: NONE/1D+5.
lineage <- data.frame(
  read_name = c("r1", "r2", "r3", "r4", "r5"),
  reference = "ref",
  ce = c("NONE_1D+5", "2D+10_1D+5", "NONE_NONE", "2D+10_NONE", "NONE_1D+5"),
  stringsAsFactors = FALSE
)

test_that("build_indel_matrix() assigns 0 to the unedited reference state", {
  im <- build_indel_matrix(lineage)
  expect_identical(im$state_map$site1[["NONE"]], 0L)
  expect_identical(im$state_map$site2[["NONE"]], 0L)
})

test_that("build_indel_matrix() assigns distinct integers to distinct indels", {
  im <- build_indel_matrix(lineage)
  # 2D+10 is a single distinct allele at site1 -> state 1, shared by r2 and r4.
  expect_identical(im$matrix["r2", "site1"], 1L)
  expect_identical(im$matrix["r2", "site1"], im$matrix["r4", "site1"])
  # r1 is reference at site1 but edited at site2.
  expect_identical(unname(im$matrix["r1", ]), c(0L, 1L))
})

test_that("build_indel_matrix() encodes missing data as -1", {
  df <- data.frame(read_name = c("a", "b"), reference = "ref",
                   ce = c("NONE_1D+5", "2D+3"), stringsAsFactors = FALSE)
  im <- build_indel_matrix(df, missing_threshold = 1)
  # b has only one site called -> site2 is missing.
  expect_identical(im$matrix["b", "site2"], -1L)
})

test_that("build_indel_matrix() drops cells above missing_threshold", {
  df <- data.frame(read_name = c("a", "b"), reference = "ref",
                   ce = c("NONE_1D+5", NA), stringsAsFactors = FALSE)
  im <- build_indel_matrix(df)  # b is all-missing -> dropped at default 0.5
  expect_true("a" %in% rownames(im$matrix))
  expect_false("b" %in% rownames(im$matrix))
})

test_that("build_indel_matrix() takes the majority call per cell", {
  df <- data.frame(
    read_name = c("x", "y", "z"),
    reference = "ref",
    ce = c("2D+10", "2D+10", "NONE"),
    stringsAsFactors = FALSE
  )
  df$cell <- "C1"  # all three reads belong to one cell
  im <- build_indel_matrix(df, cell_col = "cell")
  expect_identical(nrow(im$matrix), 1L)
  # 2D+10 appears twice vs NONE once -> majority is 2D+10 (state 1).
  expect_identical(im$matrix["C1", "site1"], 1L)
})

test_that("indel_distance() ignores missing-coded sites in pairwise comparison", {
  m <- rbind(
    a = c(0L, 1L, 2L),
    b = c(0L, 1L, -1L),   # site3 missing
    c = c(0L, 0L, 0L)
  )
  d <- as.matrix(indel_distance(m))
  # a vs b: shared sites 1,2 both equal -> distance 0.
  expect_equal(d["a", "b"], 0)
  # a vs c: sites 1 equal, 2 & 3 differ over 3 shared -> 2/3.
  expect_equal(d["a", "c"], 2 / 3)
})

test_that("indel_distance() returns 1 when cells share no non-missing site", {
  m <- rbind(a = c(0L, -1L), b = c(-1L, 0L))
  d <- as.matrix(indel_distance(m))
  expect_equal(d["a", "b"], 1)
})

test_that("indel_distance() applies per-site weights", {
  m <- rbind(a = c(0L, 0L), b = c(1L, 1L))
  # both sites differ; equal weights -> 1. Doubling site1's weight keeps it 1
  # (all differ), so test a partial-difference case instead.
  m2 <- rbind(a = c(0L, 0L), b = c(1L, 0L))  # only site1 differs
  d_equal <- as.matrix(indel_distance(m2, weights = c(1, 1)))["a", "b"]
  d_heavy <- as.matrix(indel_distance(m2, weights = c(3, 1)))["a", "b"]
  expect_equal(d_equal, 1 / 2)
  expect_equal(d_heavy, 3 / 4)
})
