.reset_cache <- function() {
  e <- asNamespace("cliqueR")$.clique_bin_cache
  e$path <- NULL
  invisible()
}

.fake_binary <- function(body = "echo \"ARGS:$@\"\nexit 0") {
  skip_on_os("windows")
  path <- tempfile(fileext = ".sh")
  writeLines(c("#!/bin/sh", body), path)
  Sys.chmod(path, mode = "0755")
  path
}

test_that(".format_clique_args handles flag types", {
  fmt <- asNamespace("cliqueR")$.format_clique_args

  expect_identical(fmt(list()), character())
  expect_identical(
    fmt(list(output_bam_file = "x.bam", threads = 4L)),
    c("--output-bam-file", "x.bam", "--threads", "4")   # snake_case -> kebab-case
  )
  expect_identical(fmt(list(correct_only = TRUE)), "--correct-only")
  expect_identical(fmt(list(correct_only = FALSE)), character())
  expect_identical(fmt(list(temp_dir = NULL)), character())
  expect_identical(fmt(list(temp_dir = NA)), character())
  expect_error(fmt(list("missing-name")), "must be named")
  expect_error(fmt(list(x = c(1, 2))), "length 1")
})

test_that("set_clique_binary validates the path", {
  .reset_cache()
  expect_error(set_clique_binary("/no/such/path"), "not found")

  non_exec <- tempfile()
  file.create(non_exec)
  Sys.chmod(non_exec, mode = "0644")
  expect_error(set_clique_binary(non_exec), "not executable")
})

test_that("clique_binary returns the path set via set_clique_binary", {
  .reset_cache()
  fake <- .fake_binary()
  set_clique_binary(fake)
  expect_identical(clique_binary(), normalizePath(fake))
})

test_that("clique_binary errors with a helpful message when nothing is findable", {
  .reset_cache()
  withr::local_envvar(CLIQUE_BIN = "")
  withr::local_dir(withr::local_tempdir())
  withr::local_path(character(), action = "replace")
  expect_error(clique_binary(refresh = TRUE), "Could not locate")
})

test_that("clique_run forwards args and captures stdout", {
  .reset_cache()
  fake <- .fake_binary()
  set_clique_binary(fake)

  res <- clique_run(
    "collapse",
    list(output_bam_file = "x.bam", correct_only = TRUE, temp_dir = NULL),
    echo = FALSE
  )
  expect_equal(res$status, 0L)
  expect_match(res$stdout, "ARGS:collapse --output-bam-file x.bam --correct-only", fixed = TRUE)
  expect_identical(res$command[1], normalizePath(fake))
})

test_that("clique_run throws on non-zero exit with stderr tail", {
  .reset_cache()
  fake <- .fake_binary("echo 'boom' >&2\nexit 7")
  set_clique_binary(fake)

  expect_error(clique_run("align", list(), echo = FALSE), "exit 7")
  expect_error(clique_run("align", list(), echo = FALSE), "boom")
})

test_that("clique_run rejects empty subcommand", {
  expect_error(clique_run("", list()), "non-empty string")
  expect_error(clique_run(c("a", "b"), list()), "single")
})
