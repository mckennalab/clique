#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L || length(arguments) > 4L) {
  stop(paste(
    "Usage: summarize-neutral-best-integration-level.R",
    "OUTPUT_ROOT LEVEL [THRESHOLD] [METHODS]"
  ))
}

root <- normalizePath(arguments[[1L]], mustWork = TRUE)
level <- suppressWarnings(as.integer(arguments[[2L]]))
threshold <- if (length(arguments) >= 3L) {
  suppressWarnings(as.numeric(arguments[[3L]]))
} else {
  0.05
}
if (is.na(level) || level < 1L || !is.finite(threshold) ||
    threshold <= 0 || threshold >= 1) {
  stop("LEVEL must be positive and THRESHOLD must be between zero and one.")
}

allowed_methods <- c(
  "nj", "parsimony", "iqtree", "mix", "vine", "cassiopeia"
)
methods <- if (length(arguments) >= 4L) {
  trimws(strsplit(arguments[[4L]], ",", fixed = TRUE)[[1L]])
} else {
  allowed_methods
}
if (!length(methods) || any(!nzchar(methods)) || anyDuplicated(methods) ||
    any(!methods %in% allowed_methods)) {
  stop("METHODS must be a unique comma-separated subset of: ",
       paste(allowed_methods, collapse = ","))
}
reconstruction_path <- file.path(
  root,
  "reconstruction",
  sprintf("level_%02d", level),
  "reconstruction_summary.csv.gz"
)
if (!file.exists(reconstruction_path)) {
  stop("Missing reconstruction summary: ", reconstruction_path)
}

results <- read.csv(gzfile(reconstruction_path), check.names = FALSE)
if (!nrow(results) || any(results$status != "success") ||
    !setequal(unique(results$method), methods) ||
    anyDuplicated(results[c("condition_id", "method")])) {
  stop("The reconstruction level is incomplete or contains failed/duplicate results.")
}
effective_level <- ifelse(
  results$system == "mitochondrial",
  results$observation_depth,
  results$integrations
)
if (anyNA(effective_level) || any(effective_level != level)) {
  stop("The reconstruction summary contains a condition from another level.")
}
if (any(results$shape != "neutral") || length(unique(results$tree_size)) != 1L) {
  stop("Only one neutral-tree sample size may be summarized at a time.")
}

confidence_summary <- function(values_by_seed) {
  replicate_means <- vapply(values_by_seed, mean, numeric(1))
  replicate_count <- length(replicate_means)
  replicate_sd <- if (replicate_count > 1L) {
    stats::sd(replicate_means)
  } else {
    NA_real_
  }
  replicate_se <- replicate_sd / sqrt(replicate_count)
  radius <- if (replicate_count > 1L) {
    stats::qt(0.975, df = replicate_count - 1L) * replicate_se
  } else {
    NA_real_
  }
  mean_value <- mean(replicate_means)
  c(
    replicates = replicate_count,
    mean_normalized_rf = mean_value,
    replicate_sd_normalized_rf = replicate_sd,
    replicate_se_normalized_rf = replicate_se,
    ci95_lower_normalized_rf = if (is.na(radius)) NA_real_ else {
      max(0, mean_value - radius)
    },
    ci95_upper_normalized_rf = if (is.na(radius)) NA_real_ else {
      min(1, mean_value + radius)
    }
  )
}

summarize_part <- function(part) {
  confidence_summary(split(part$normalized_rf, part$population_seed))
}

method_rows <- list()
method_index <- 0L
for (system in sort(unique(results$system))) {
  system_results <- results[results$system == system, , drop = FALSE]
  for (method in methods) {
    part <- system_results[system_results$method == method, , drop = FALSE]
    summary_values <- summarize_part(part)
    method_index <- method_index + 1L
    method_rows[[method_index]] <- data.frame(
      system = system,
      level = level,
      method = method,
      trees = nrow(part),
      as.list(summary_values),
      stringsAsFactors = FALSE
    )
  }
}
method_summary <- do.call(rbind, method_rows)
rownames(method_summary) <- NULL

overall_rows <- lapply(sort(unique(results$system)), function(system) {
  part <- results[results$system == system, , drop = FALSE]
  summary_values <- summarize_part(part)
  method_part <- method_summary[method_summary$system == system, , drop = FALSE]
  best <- method_part[which.min(method_part$mean_normalized_rf), , drop = FALSE]
  data.frame(
    system = system,
    level = level,
    tree_size = unique(part$tree_size),
    trees = nrow(part),
    as.list(summary_values),
    best_method = best$method,
    best_method_mean_normalized_rf = best$mean_normalized_rf,
    threshold = threshold,
    threshold_reached = summary_values[["mean_normalized_rf"]] < threshold,
    stringsAsFactors = FALSE
  )
})
overall_summary <- do.call(rbind, overall_rows)
rownames(overall_summary) <- NULL

replace_level_rows <- function(path, current) {
  if (file.exists(path)) {
    previous <- read.delim(path, check.names = FALSE)
    current_keys <- paste(current$system, current$level, sep = "\r")
    previous_keys <- paste(previous$system, previous$level, sep = "\r")
    previous <- previous[!previous_keys %in% current_keys, , drop = FALSE]
    current <- rbind(previous[, names(current), drop = FALSE], current)
  }
  current <- current[order(current$system, current$level), , drop = FALSE]
  write.table(current, path, sep = "\t", row.names = FALSE, quote = FALSE)
  current
}

overall_path <- file.path(root, "accuracy_by_recorder_integration.tsv")
method_path <- file.path(root, "accuracy_by_recorder_integration_method.tsv")
all_overall <- replace_level_rows(overall_path, overall_summary)
invisible(replace_level_rows(method_path, method_summary))

current_active <- overall_summary$system[!overall_summary$threshold_reached]
writeLines(
  paste(current_active, collapse = ","),
  file.path(root, "active_systems.txt")
)

cat(sprintf(
  "Level %d: %d recorder(s) tested; %d crossed normalized RF < %.3f.\n",
  level,
  nrow(overall_summary),
  sum(overall_summary$threshold_reached),
  threshold
))
print(overall_summary[c(
  "system", "level", "mean_normalized_rf", "best_method",
  "best_method_mean_normalized_rf", "threshold_reached"
)], row.names = FALSE)
cat("ACTIVE_SYSTEMS=", paste(current_active, collapse = ","), "\n", sep = "")
