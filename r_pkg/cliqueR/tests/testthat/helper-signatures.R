expect_args <- function(fn, expected) {
  expect_true(is.function(fn))
  expect_identical(names(formals(fn)), expected)
}
