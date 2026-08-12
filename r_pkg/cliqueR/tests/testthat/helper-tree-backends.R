make_lineage_matrix <- function() {
  m <- rbind(
    "cell:A" = c(0L, 0L, 0L),
    "cell B" = c(1L, 0L, 1L),
    "cell/C" = c(2L, 1L, -1L),
    "cell_D" = c(1L, 1L, 0L)
  )
  colnames(m) <- c("target_1", "target_2", "target_3")
  m
}

make_fake_executable <- function(lines) {
  path <- tempfile("cliqueR-fake-")
  writeLines(c("#!/bin/sh", "set -eu", lines), path, useBytes = TRUE)
  Sys.chmod(path, mode = "0755")
  path
}

safe_four_tip_newick <- function() {
  "((clq0000001:1,clq0000002:1):1,(clq0000003:1,clq0000004:1):1);"
}

cassiopeia_test_available <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  python <- unname(Sys.which(c("python3", "python")))
  python <- python[nzchar(python)][1L]
  if (length(python) == 0L || is.na(python)) return(FALSE)
  withr::local_envvar(RETICULATE_PYTHON = python)
  tryCatch({
    suppressWarnings(reticulate::import("cassiopeia", delay_load = FALSE))
    TRUE
  }, error = function(e) FALSE)
}
