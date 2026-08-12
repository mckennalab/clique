.format_clique_args <- function(args) {
  if (length(args) == 0) return(character())
  nms <- names(args)
  if (is.null(nms) || any(!nzchar(nms))) {
    rlang::abort("All elements of `args` must be named.")
  }
  out <- character()
  for (nm in nms) {
    val <- args[[nm]]
    if (is.null(val)) next
    if (length(val) == 1 && is.na(val)) next
    # The clique CLI (clap) exposes kebab-case flags derived from its field
    # names (e.g. `--read-structure`), so translate the idiomatic snake_case
    # names R callers pass. Names given with dashes already pass through.
    flag <- paste0("--", gsub("_", "-", nm))
    if (isTRUE(val)) {
      out <- c(out, flag)
    } else if (isFALSE(val)) {
      next
    } else {
      if (length(val) != 1) {
        rlang::abort(sprintf("Argument `%s` must be length 1, got %d.", nm, length(val)))
      }
      out <- c(out, flag, as.character(val))
    }
  }
  out
}

#' Low-level invocation of the clique CLI
#'
#' Runs a single `clique` subcommand and captures stdout/stderr. Used internally
#' by [align()] and [collapse()]; exported for users who need to call
#' subcommands or flags this package doesn't yet wrap.
#'
#' Argument conventions for `args`:
#' * `TRUE` becomes a bare `--flag`.
#' * `FALSE`, `NULL`, and `NA` are dropped (so callers can pass optional flags
#'   without conditional logic).
#' * Other scalar values become `--flag value` with `as.character()` coercion.
#' * Element names are translated to the CLI's kebab-case flags: underscores
#'   become dashes, so `output_bam_file` is sent as `--output-bam-file`.
#'
#' @param subcommand Name of the clique subcommand (e.g. `"align"`, `"collapse"`).
#' @param args Named list of `--flag = value` pairs (see details).
#' @param echo If `TRUE`, stream stdout/stderr to the R console live.
#' @param timeout Timeout in seconds (`Inf` for none).
#' @return Invisibly, a list with `status`, `stdout`, `stderr`, and `command`
#'   (the resolved argv).
#' @examples
#' \dontrun{
#' # Equivalent to: clique collapse --input_bam_file in.bam \
#' #   --output_bam_file out.bam --read_structure rs.yaml --correct_only
#' clique_run("collapse", list(
#'   input_bam_file = "in.bam",
#'   output_bam_file = "out.bam",
#'   read_structure = "rs.yaml",
#'   correct_only = TRUE
#' ))
#'
#' # Print the version banner
#' clique_run("--version", list())
#' }
#' @export
clique_run <- function(subcommand, args = list(), echo = TRUE, timeout = Inf) {
  if (!is.character(subcommand) || length(subcommand) != 1 || !nzchar(subcommand)) {
    rlang::abort("`subcommand` must be a single non-empty string.")
  }
  bin <- clique_binary()
  argv <- c(subcommand, .format_clique_args(args))

  result <- processx::run(
    command = bin,
    args = argv,
    echo = echo,
    echo_cmd = FALSE,
    spinner = FALSE,
    error_on_status = FALSE,
    timeout = timeout
  )

  out <- list(
    status = result$status,
    stdout = result$stdout,
    stderr = result$stderr,
    command = c(bin, argv)
  )

  if (!identical(result$status, 0L)) {
    stderr_lines <- strsplit(result$stderr %||% "", "\n", fixed = TRUE)[[1]]
    tail_lines <- utils::tail(stderr_lines, 20)
    rlang::abort(c(
      sprintf("clique %s failed (exit %s)", subcommand, format(result$status)),
      i = paste(tail_lines, collapse = "\n")
    ))
  }

  invisible(out)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
