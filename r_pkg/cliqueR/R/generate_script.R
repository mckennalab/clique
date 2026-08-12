#' Require that a file exists, with a clear error.
#' @noRd
.require_file <- function(path, what) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    rlang::abort(sprintf("%s not found: %s", what, if (is.null(path)) "<NULL>" else path))
  }
  invisible(path)
}

#' Extract the reference names (keys) from a clique read-structure YAML without
#' a YAML parser (avoids choking on the `!Read1`-style custom tags).
#' @noRd
.yaml_reference_names <- function(yaml_path) {
  lines <- readLines(yaml_path, warn = FALSE)
  ref_start <- grep("^references:\\s*$", lines)
  if (length(ref_start) == 0) return(character())
  names <- character()
  tail_lines <- lines[(ref_start[1] + 1):length(lines)]
  for (ln in tail_lines) {
    if (grepl("^\\S", ln)) break                       # dedented to column 0 -> block ended
    if (grepl("^  [^ ].*:\\s*$", ln)) {                # exactly-2-space-indented `name:`
      key <- sub(":\\s*$", "", sub("^  ", "", ln))
      names <- c(names, gsub('^"|"$', "", trimws(key)))
    }
  }
  names
}

#' Warn if the FASTA and YAML disagree on reference names.
#' @noRd
.check_reference_names <- function(reference, read_structure) {
  fasta_names <- .read_fasta(reference)$name
  yaml_names <- .yaml_reference_names(read_structure)
  if (length(yaml_names) == 0) {
    rlang::warn("Could not find a `references:` block in the read-structure YAML to cross-check.")
    return(invisible())
  }
  only_fasta <- setdiff(fasta_names, yaml_names)
  only_yaml <- setdiff(yaml_names, fasta_names)
  if (length(only_fasta) || length(only_yaml)) {
    msg <- "Reference names differ between the FASTA and the read-structure YAML:"
    if (length(only_fasta)) msg <- c(msg, i = sprintf("in FASTA only: %s", paste(only_fasta, collapse = ", ")))
    if (length(only_yaml)) msg <- c(msg, i = sprintf("in YAML only: %s", paste(only_yaml, collapse = ", ")))
    rlang::warn(msg)
  }
  invisible()
}

#' Generate a bash script for the clique align (+ collapse) pipeline
#'
#' Writes a self-contained bash script that runs `clique align` on the given read
#' structure and reads -- aligning to the references embedded in the read
#' structure and calling edit events (the `ce` BAM tag) -- then sorts and indexes
#' the BAM, and optionally runs `clique collapse` for per-molecule consensus.
#'
#' The reference FASTA is not a direct `clique` input (references live inside the
#' read-structure YAML). It is used here to sanity-check that the YAML's
#' reference names match the FASTA, and is recorded in the script header.
#'
#' @param read_structure Path to the read-structure YAML.
#' @param reference Path to the reference FASTA (for the name cross-check).
#' @param read1 Path to the read-1 FASTQ (required).
#' @param output_script Path to write the bash script.
#' @param read2,index1,index2 Optional additional read/index FASTQs.
#' @param output_dir Directory the script writes BAMs into (it is created).
#' @param sample_name Base name for outputs; defaults to the `read1` basename.
#' @param clique_bin Path to the `clique` executable; defaults to
#'   [clique_binary()] when resolvable, else `"clique"` on `PATH`. The script also
#'   honours a `CLIQUE_BIN` environment variable at run time.
#' @param samtools Path to `samtools`.
#' @param threads Threads for align / collapse / sort.
#' @param aligner clique `--aligner` value (`"WFA"`, `"Degenerate"`, `"Inversion"`).
#' @param min_read_length,max_reference_multiplier clique align tuning flags.
#' @param collapse If `TRUE`, also emit a `clique collapse` step (forces
#'   `sort_and_index`).
#' @param correct_only Passed to collapse as `--correct-only` (correct tags only,
#'   no consensus).
#' @param sort_and_index Sort + index the aligned BAM (required before collapse,
#'   and for IGV).
#' @param validate_reference Cross-check YAML reference names against the FASTA.
#' @param execute If `TRUE`, run the generated script with `bash` via [processx].
#' @return Invisibly, the path to the written (and, if `execute`, run) script.
#' @details
#' With `sample_name = "sample"`, the script writes `sample.aligned.bam`,
#' optionally `sample.sorted.bam` plus its index, and, when `collapse = TRUE`,
#' `sample.consensus.bam` under `output_dir`. The generated script requires
#' `bash`, the `clique` binary, and `samtools` when sorting is enabled.
#' @seealso [fasta_targets_to_yaml()], [clique_run()]
#' @examples
#' \dontrun{
#' generate_clique_script(
#'   read_structure = "read_structure.yaml",
#'   reference      = "amplicons.fa",
#'   read1          = "sample_R1.fastq.gz",
#'   read2          = "sample_R2.fastq.gz",
#'   output_script  = "run_clique.sh",
#'   collapse       = TRUE
#' )
#' }
#' @export
generate_clique_script <- function(read_structure,
                                   reference,
                                   read1,
                                   output_script,
                                   read2 = NULL,
                                   index1 = NULL,
                                   index2 = NULL,
                                   output_dir = ".",
                                   sample_name = NULL,
                                   clique_bin = NULL,
                                   samtools = "samtools",
                                   threads = 4L,
                                   aligner = "WFA",
                                   min_read_length = 50L,
                                   max_reference_multiplier = 2L,
                                   collapse = FALSE,
                                   correct_only = FALSE,
                                   sort_and_index = TRUE,
                                   validate_reference = TRUE,
                                   execute = FALSE) {

  .require_file(read_structure, "read_structure (YAML)")
  .require_file(reference, "reference (FASTA)")
  .require_file(read1, "read1 FASTQ")
  if (!is.null(read2)) .require_file(read2, "read2 FASTQ")
  if (!is.null(index1)) .require_file(index1, "index1 FASTQ")
  if (!is.null(index2)) .require_file(index2, "index2 FASTQ")

  if (!tolower(aligner) %in% c("wfa", "degenerate", "inversion")) {
    rlang::abort(sprintf("Unknown aligner '%s' (one of WFA, Degenerate, Inversion).", aligner))
  }
  if (isTRUE(collapse) && !isTRUE(sort_and_index)) {
    rlang::inform("collapse = TRUE requires a sorted, indexed BAM; enabling sort_and_index.")
    sort_and_index <- TRUE
  }

  bin <- clique_bin
  if (is.null(bin)) bin <- tryCatch(clique_binary(), error = function(e) NULL)
  if (is.null(bin)) {
    bin <- "clique"
    rlang::warn("Could not resolve the clique binary; the script will use `clique` from PATH (override with CLIQUE_BIN or clique_bin=).")
  }

  if (isTRUE(validate_reference)) .check_reference_names(reference, read_structure)

  if (is.null(sample_name)) {
    sample_name <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(read1))
    sample_name <- sub("[._]R?1(_001)?$", "", sample_name)
  }

  abspath <- function(p) if (is.null(p)) NULL else normalizePath(p, mustWork = FALSE)
  yaml_p  <- abspath(read_structure)
  ref_p   <- abspath(reference)
  r1_p    <- abspath(read1)
  r2_p    <- abspath(read2)
  i1_p    <- abspath(index1)
  i2_p    <- abspath(index2)
  out_p   <- abspath(output_dir)

  aligned <- sprintf('"$OUTDIR/%s.aligned.bam"', sample_name)
  sorted  <- sprintf('"$OUTDIR/%s.sorted.bam"', sample_name)
  consensus <- sprintf('"$OUTDIR/%s.consensus.bam"', sample_name)

  # --- assemble the align command (optional flags only when supplied) ---
  align <- c(
    '"$CLIQUE" align \\',
    sprintf('  --read-structure "%s" \\', yaml_p),
    sprintf('  --read1 "%s" \\', r1_p)
  )
  if (!is.null(r2_p)) align <- c(align, sprintf('  --read2 "%s" \\', r2_p))
  if (!is.null(i1_p)) align <- c(align, sprintf('  --index1 "%s" \\', i1_p))
  if (!is.null(i2_p)) align <- c(align, sprintf('  --index2 "%s" \\', i2_p))
  align <- c(align,
    sprintf('  --output-bam-file %s \\', aligned),
    sprintf('  --aligner %s \\', tolower(aligner)),
    sprintf('  --threads %d \\', as.integer(threads)),
    sprintf('  --min-read-length %d \\', as.integer(min_read_length)),
    sprintf('  --max-reference-multiplier %d', as.integer(max_reference_multiplier))
  )

  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  script <- c(
    "#!/usr/bin/env bash",
    sprintf("# Generated by cliqueR::generate_clique_script (%s)", ts),
    sprintf("# read structure : %s", yaml_p),
    sprintf("# reference       : %s", ref_p),
    "# Aligns reads to the references embedded in the read structure and calls",
    "# edit events (the `ce` BAM tag); then sorts/indexes; then optionally collapses.",
    "set -euo pipefail",
    "",
    sprintf('CLIQUE="${CLIQUE_BIN:-%s}"', bin),
    sprintf('SAMTOOLS="%s"', samtools),
    sprintf('OUTDIR="%s"', out_p),
    'mkdir -p "$OUTDIR"',
    "",
    "echo '[clique] 1/3 aligning reads and calling events...'",
    align,
    ""
  )

  if (isTRUE(sort_and_index)) {
    script <- c(script,
      "echo '[clique] 2/3 sorting and indexing...'",
      sprintf('"$SAMTOOLS" sort -@ %d -o %s %s', as.integer(threads), sorted, aligned),
      sprintf('"$SAMTOOLS" index %s', sorted),
      ""
    )
  }

  if (isTRUE(collapse)) {
    collapse_cmd <- c(
      '"$CLIQUE" collapse \\',
      sprintf('  --read-structure "%s" \\', yaml_p),
      sprintf('  --input-bam-file %s \\', sorted),
      sprintf('  --output-bam-file %s \\', consensus),
      sprintf('  --threads %d%s', as.integer(threads), if (isTRUE(correct_only)) " \\" else "")
    )
    if (isTRUE(correct_only)) collapse_cmd <- c(collapse_cmd, "  --correct-only")
    script <- c(script,
      "echo '[clique] 3/3 collapsing by UMI to per-molecule consensus...'",
      collapse_cmd,
      ""
    )
  }

  final_bam <- if (isTRUE(collapse)) consensus else if (isTRUE(sort_and_index)) sorted else aligned
  final_bam_disp <- gsub('^"|"$', "", final_bam)   # strip the outer quotes for display
  script <- c(script, sprintf('echo "[clique] done -> %s"', final_bam_disp))

  writeLines(script, output_script)
  Sys.chmod(output_script, mode = "0755")
  rlang::inform(sprintf("Wrote %s (bash). Run it with: bash %s", output_script, output_script))

  if (isTRUE(execute)) {
    rlang::inform("Executing the generated script...")
    processx::run("bash", output_script, echo = TRUE, error_on_status = TRUE)
  }

  invisible(output_script)
}
