.clique_bin_cache <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .clique_bin_cache$path <- NULL
  invisible()
}

.clique_exe_name <- function() {
  if (.Platform$OS.type == "windows") "clique.exe" else "clique"
}

.is_executable <- function(path) {
  file.exists(path) && file.access(path, mode = 1L) == 0L
}

#' Locate the clique binary
#'
#' Resolution order:
#' 1. `CLIQUE_BIN` environment variable, if set and executable.
#' 2. `clique` on the user's `PATH` (via [Sys.which()]).
#' 3. `rust_cmd/target/release/clique` relative to the current working
#'    directory, for developers building from this monorepo.
#'
#' The resolved path is cached for the session. Call with `refresh = TRUE`
#' to re-resolve (e.g. after rebuilding the binary).
#'
#' @param refresh If `TRUE`, ignore the cached path and re-resolve.
#' @return Absolute path to the `clique` executable.
#' @examples
#' \dontrun{
#' # Default: looks at CLIQUE_BIN, then PATH, then rust_cmd/target/release/
#' clique_binary()
#'
#' # Force re-resolution after rebuilding the Rust binary
#' clique_binary(refresh = TRUE)
#' }
#' @export
clique_binary <- function(refresh = FALSE) {
  if (!refresh && !is.null(.clique_bin_cache$path)) {
    return(.clique_bin_cache$path)
  }

  exe <- .clique_exe_name()
  candidates <- character()

  env_bin <- Sys.getenv("CLIQUE_BIN", unset = "")
  if (nzchar(env_bin)) candidates <- c(candidates, env_bin)

  on_path <- unname(Sys.which(exe))
  if (nzchar(on_path)) candidates <- c(candidates, on_path)

  local_build <- file.path("rust_cmd", "target", "release", exe)
  if (file.exists(local_build)) candidates <- c(candidates, local_build)

  for (cand in candidates) {
    if (.is_executable(cand)) {
      .clique_bin_cache$path <- normalizePath(cand, mustWork = TRUE)
      return(.clique_bin_cache$path)
    }
  }

  rlang::abort(c(
    "Could not locate the `clique` binary.",
    i = "Set the CLIQUE_BIN environment variable, or call `set_clique_binary(path)`.",
    i = "To build it: cd rust_cmd && cargo build --release",
    i = sprintf("Checked: %s", paste(candidates, collapse = ", "))
  ))
}

#' Override the clique binary location for the current session
#'
#' Sets the cached binary path used by [clique_binary()]. The path must exist
#' and be executable.
#'
#' @param path Path to a `clique` executable.
#' @return Invisibly, the normalized path that was set.
#' @examples
#' \dontrun{
#' set_clique_binary("~/code/clique/rust_cmd/target/release/clique")
#' clique_binary()  # confirms the override
#' }
#' @export
set_clique_binary <- function(path) {
  if (!file.exists(path)) {
    rlang::abort(sprintf("Binary not found: %s", path))
  }
  if (!.is_executable(path)) {
    rlang::abort(sprintf("File is not executable: %s", path))
  }
  .clique_bin_cache$path <- normalizePath(path, mustWork = TRUE)
  invisible(.clique_bin_cache$path)
}
