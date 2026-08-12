# clique target-type names accepted by the read-structure YAML.
.clique_target_types <- c(
  "Static", "Cas9WT", "Cas12AWT", "Cas9ABE", "Cas9CBE", "Cas9ABECBE",
  "Cas12ABE", "Cas12CBE", "Cas12ABECBE", "Cas9Homing", "Cas9ABEPalindrome"
)
.clique_orientations <- c("Forward", "Reverse", "ReverseComplement", "Unknown")
.clique_read_types <- c("Read1", "Read2", "Index1", "Index2", "Spacer")
.clique_merge_strategies <- c("Align", "Concatenate", "ConcatenateBothForward")

#' Read a (multi-record) FASTA into a data frame of name/sequence.
#' @noRd
.read_fasta <- function(path) {
  if (!file.exists(path)) rlang::abort(sprintf("FASTA file not found: %s", path))
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*$", lines)]                 # drop blank lines
  if (length(lines) == 0) rlang::abort("FASTA file is empty.")
  hdr <- grep("^>", lines)
  if (length(hdr) == 0) rlang::abort("No FASTA headers ('>') found.")

  names <- sub("\\s.*$", "", sub("^>\\s*", "", lines[hdr]))  # first token after '>'
  ends <- c(hdr[-1] - 1L, length(lines))
  seqs <- character(length(hdr))
  for (i in seq_along(hdr)) {
    body <- if (hdr[i] + 1L <= ends[i]) lines[(hdr[i] + 1L):ends[i]] else character()
    seqs[i] <- toupper(gsub("[^A-Za-z]", "", paste0(body, collapse = "")))
  }

  if (any(!nzchar(names))) rlang::abort("A FASTA record has an empty name.")
  if (any(!nzchar(seqs))) {
    rlang::abort(sprintf("FASTA record(s) with empty sequence: %s",
                         paste(names[!nzchar(seqs)], collapse = ", ")))
  }
  if (anyDuplicated(names)) {
    dups <- unique(names[duplicated(names)])
    rlang::abort(sprintf("Duplicate reference name(s) in FASTA: %s", paste(dups, collapse = ", ")))
  }
  data.frame(name = names, sequence = seqs, stringsAsFactors = FALSE)
}

#' Reverse complement an (upper-cased) IUPAC DNA string.
#' @noRd
.revcomp <- function(s) {
  comp <- chartr("ACGTUNRYSWKMBDHV", "TGCAANYRSWMKVHDB", toupper(s))
  paste(rev(strsplit(comp, "", fixed = TRUE)[[1]]), collapse = "")
}

#' Normalize the `targets` argument into a data frame with `target` and `type`.
#' @noRd
.normalize_targets <- function(targets, target_types) {
  if (is.data.frame(targets)) {
    tcol <- intersect(c("target", "sequence", "seq"), names(targets))[1]
    ycol <- intersect(c("type", "target_type"), names(targets))[1]
    if (is.na(tcol)) rlang::abort("data.frame `targets` must have a `target` (or `sequence`) column.")
    type <- if (!is.na(ycol)) as.character(targets[[ycol]]) else rep(target_types, length.out = nrow(targets))
    out <- data.frame(target = as.character(targets[[tcol]]), type = type, stringsAsFactors = FALSE)
  } else {
    targets <- as.character(targets)
    if (length(target_types) == 1L) target_types <- rep(target_types, length(targets))
    if (length(target_types) != length(targets)) {
      rlang::abort(sprintf("`target_types` length %d must be 1 or equal to `targets` length %d.",
                           length(target_types), length(targets)))
    }
    out <- data.frame(target = targets, type = as.character(target_types), stringsAsFactors = FALSE)
  }
  if (nrow(out) == 0L) rlang::abort("No targets supplied.")
  if (any(!nzchar(out$target)) || any(!grepl("^[A-Za-z]+$", out$target))) {
    rlang::abort("Every target must be a non-empty string of DNA letters.")
  }
  out
}

#' Validate the `reads` layout list.
#' @noRd
.validate_reads <- function(reads) {
  if (!is.list(reads) || length(reads) == 0L) rlang::abort("`reads` must be a non-empty list.")
  for (r in reads) {
    if (is.null(r$type) || !r$type %in% .clique_read_types) {
      rlang::abort(sprintf("Each read needs a valid `type` (one of %s).",
                           paste(.clique_read_types, collapse = ", ")))
    }
    if (identical(r$type, "Spacer")) {
      if (is.null(r$spacer_sequence) || !nzchar(r$spacer_sequence)) {
        rlang::abort("A Spacer read needs a non-empty `spacer_sequence`.")
      }
    } else {
      o <- if (is.null(r$orientation)) "Forward" else r$orientation
      if (!o %in% .clique_orientations) {
        rlang::abort(sprintf("Invalid read orientation '%s' (one of %s).",
                             o, paste(.clique_orientations, collapse = ", ")))
      }
    }
  }
}

#' Quote a YAML key only if it isn't a plain identifier.
#' @noRd
.yaml_key <- function(k) if (grepl("^[A-Za-z0-9_.-]+$", k)) k else paste0('"', gsub('"', '\\\\"', k), '"')

#' Render the read-structure YAML text.
#' @noRd
.render_read_structure_yaml <- function(refs, matched, merge, known_strand, reads, uppercase) {
  q <- function(s) paste0('"', s, '"')
  arr <- function(v) if (length(v) == 0L) "[]" else paste0("[", paste(vapply(v, q, character(1)), collapse = ", "), "]")

  lines <- "---"
  if (!is.null(merge)) lines <- c(lines, sprintf("merge: %s", q(merge)))
  lines <- c(lines, sprintf("known_strand: %s", if (isTRUE(known_strand)) "true" else "false"), "reads:")
  for (r in reads) {
    if (identical(r$type, "Spacer")) {
      lines <- c(lines, "  - !Spacer", sprintf("    spacer_sequence: %s", q(r$spacer_sequence)))
    } else {
      o <- if (is.null(r$orientation)) "Forward" else r$orientation
      lines <- c(lines, sprintf("  - !%s", r$type), sprintf("    orientation: %s", o))
    }
  }
  lines <- c(lines, "references:")
  for (i in seq_len(nrow(refs))) {
    seqi <- if (uppercase) toupper(refs$sequence[i]) else refs$sequence[i]
    m <- matched[[i]]
    tgt <- if (uppercase) toupper(m$target) else m$target
    lines <- c(lines,
      sprintf("  %s:", .yaml_key(refs$name[i])),
      sprintf("    sequence: %s", q(seqi)),
      sprintf("    targets: %s", arr(tgt)),
      sprintf("    target_types: %s", arr(m$type)),
      "    umi_configurations: {}")
  }
  paste(lines, collapse = "\n")
}

#' Build a clique read-structure YAML from a multi-reference FASTA + target set
#'
#' Reads a (multi-record) FASTA and a set of CRISPR target/protospacer sequences,
#' locates each target within every reference, and writes a read-structure YAML
#' that `clique align` / `clique collapse` accept. Every target written for a
#' reference is guaranteed to be a forward-strand substring of that reference's
#' sequence (clique's loader requires this), so targets that only occur on the
#' reverse strand are written as their reverse complement.
#'
#' @param fasta_file Path to a FASTA file with one or more reference records.
#' @param targets Target sequences: a character vector, or a data frame with a
#'   `target` (or `sequence`) column and optionally a `type` column.
#' @param output_file Path to write the YAML.
#' @param target_types clique target type(s) used when `targets` is a plain
#'   character vector (or a data frame without a `type` column): a single value
#'   recycled to all, or a vector the same length as `targets`. Must be one of
#'   the clique target types (e.g. `"Cas9WT"`, `"Cas9ABE"`, `"Cas9ABEPalindrome"`).
#'   `PrimeEdit` is not accepted by this helper because it also requires a
#'   `prime_edits` specification; add prime-edit targets directly to the YAML.
#' @param search_reverse_complement If `TRUE`, a target not found on the forward
#'   strand is searched as its reverse complement; when found that way, the
#'   reverse complement (the actual forward-strand substring) is written.
#' @param keep_references_without_targets Keep references in which no target
#'   matched (emitting empty `targets`). If `FALSE`, such references are dropped.
#' @param merge clique merge strategy (`"Align"`, `"Concatenate"`,
#'   `"ConcatenateBothForward"`) or `NULL` to omit.
#' @param known_strand Value of the layout's `known_strand` flag.
#' @param reads A list describing read positions; each element is a list with
#'   `type` (`Read1`/`Read2`/`Index1`/`Index2`/`Spacer`) plus, for reads,
#'   `orientation` (default `"Forward"`) or, for a spacer, `spacer_sequence`.
#' @param uppercase Upper-case reference and target sequences in the output
#'   (recommended; clique matches targets case-sensitively).
#' @return Invisibly, a named list (per kept reference) of the matched targets
#'   with their type and strand (`"+"`/`"-"`). The YAML is written to
#'   `output_file`.
#' @seealso [generate_clique_script()]
#' @examples
#' \dontrun{
#' fasta_targets_to_yaml(
#'   fasta_file   = "amplicons.fa",
#'   targets      = c("GACGGCTATACAAGGCATCGCGG", "CTCGTCAATACACCTTACGGAGG"),
#'   output_file  = "read_structure.yaml",
#'   target_types = "Cas9ABE"
#' )
#' }
#' @export
fasta_targets_to_yaml <- function(fasta_file,
                                  targets,
                                  output_file,
                                  target_types = "Cas9WT",
                                  search_reverse_complement = TRUE,
                                  keep_references_without_targets = TRUE,
                                  merge = "ConcatenateBothForward",
                                  known_strand = TRUE,
                                  reads = list(list(type = "Read1", orientation = "Forward")),
                                  uppercase = TRUE) {

  tdf <- .normalize_targets(targets, target_types)
  bad_types <- setdiff(unique(tdf$type), .clique_target_types)
  if (length(bad_types)) {
    rlang::abort(c(sprintf("Unknown clique target type(s): %s", paste(bad_types, collapse = ", ")),
                  i = paste("Valid types:", paste(.clique_target_types, collapse = ", "))))
  }
  tdf$target <- toupper(tdf$target)

  if (!is.null(merge) && !merge %in% .clique_merge_strategies) {
    rlang::abort(sprintf("Unknown merge strategy '%s' (one of %s).",
                         merge, paste(.clique_merge_strategies, collapse = ", ")))
  }
  .validate_reads(reads)

  refs <- .read_fasta(fasta_file)

  matched <- vector("list", nrow(refs))
  found_any <- rep(FALSE, nrow(tdf))
  for (i in seq_len(nrow(refs))) {
    seqi <- refs$sequence[i]
    hit_target <- character(0); hit_type <- character(0); hit_strand <- character(0)
    for (j in seq_len(nrow(tdf))) {
      tj <- tdf$target[j]
      if (grepl(tj, seqi, fixed = TRUE)) {
        hit_target <- c(hit_target, tj); hit_type <- c(hit_type, tdf$type[j]); hit_strand <- c(hit_strand, "+")
        found_any[j] <- TRUE
      } else if (isTRUE(search_reverse_complement)) {
        rc <- .revcomp(tj)
        if (grepl(rc, seqi, fixed = TRUE)) {
          hit_target <- c(hit_target, rc); hit_type <- c(hit_type, tdf$type[j]); hit_strand <- c(hit_strand, "-")
          found_any[j] <- TRUE
        }
      }
    }
    keep <- !duplicated(hit_target)   # clique locates each target once (first occurrence)
    matched[[i]] <- list(target = hit_target[keep], type = hit_type[keep], strand = hit_strand[keep])
  }

  if (any(!found_any)) {
    nf <- tdf$target[!found_any]
    message(sprintf("%d of %d target(s) matched no reference%s", sum(!found_any), nrow(tdf),
                    if (length(nf) <= 8) paste0(": ", paste(nf, collapse = ", ")) else ""))
  }

  keep_ref <- if (isTRUE(keep_references_without_targets)) {
    rep(TRUE, nrow(refs))
  } else {
    vapply(matched, function(m) length(m$target) > 0L, logical(1))
  }
  if (!any(keep_ref)) {
    rlang::abort("No references retained (no target matched and keep_references_without_targets = FALSE).")
  }

  yaml <- .render_read_structure_yaml(
    refs[keep_ref, , drop = FALSE], matched[keep_ref], merge, known_strand, reads, uppercase
  )
  writeLines(yaml, output_file)

  n_per_ref <- vapply(matched[keep_ref], function(m) length(m$target), integer(1))
  message(sprintf("Wrote %s: %d reference(s), %d target placement(s).",
                  output_file, sum(keep_ref), sum(n_per_ref)))

  invisible(stats::setNames(matched[keep_ref], refs$name[keep_ref]))
}
