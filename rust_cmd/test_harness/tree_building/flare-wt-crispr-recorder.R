# FLARE-like wild-type Cas9 recorder support for the fixed-tree benchmark.
#
# This extension uses remote_mito's ordinary nuclease event model and adds the
# inter-target deletion behavior needed by a multi-cut WT-CRISPR array. It is
# intentionally a benchmark model, not a claim of fitted FLARE parameters.

flare_wt_crispr_probability <- function(value, name, allow_one = FALSE) {
  value <- suppressWarnings(as.numeric(value))
  invalid_upper <- if (allow_one) value > 1 else value >= 1
  if (length(value) != 1L || !is.finite(value) || value < 0 || invalid_upper) {
    stop(
      name,
      if (allow_one) " must be one probability in [0, 1]." else
        " must be one probability in [0, 1)."
    )
  }
  value
}

configure_flare_wt_crispr_params <- function(
    params,
    insertion_probability,
    deletion_probability,
    dropout_probability = 0.30,
    dropout_radius = 200L,
    deletion_extent = NULL) {
  insertion_probability <- flare_wt_crispr_probability(
    insertion_probability,
    "WT-CRISPR insertion probability"
  )
  deletion_probability <- flare_wt_crispr_probability(
    deletion_probability,
    "WT-CRISPR deletion probability"
  )
  dropout_probability <- flare_wt_crispr_probability(
    dropout_probability,
    "WT-CRISPR inter-target dropout probability",
    allow_one = TRUE
  )
  dropout_radius <- suppressWarnings(as.integer(dropout_radius))
  if (length(dropout_radius) != 1L || is.na(dropout_radius) ||
      dropout_radius < 0L) {
    stop("WT-CRISPR inter-target dropout radius must be non-negative.")
  }

  # Ten targets at 15, 36, ..., 204 reproduce the documented S:15:20 array.
  params$bc_length <- 215L
  params$physicell_adapter$recorder_system <- "FLARE-like WT CRISPR/Cas9"
  params$be_targets <- lineage_benchmark_empty_target_spec()
  params$nuclease_targets <- list(
    num_targets = 10L,
    edit_rate_class_fractions = list(high = 1, medium = 0, low = 0),
    config = "S:15:20",
    editing_window = list(
      size = 0L, decaying = FALSE, close_after_edit = TRUE
    ),
    interdeletion_dropout_probability = dropout_probability,
    interdeletion_dropout_radius = dropout_radius,
    prime_editing_system = FALSE
  )
  params$prime_editing_backend <- NULL
  params$palincode_adapter <- NULL
  params$flare_wt_crispr_backend <- list(
    enabled = TRUE,
    insertion_probability_per_target_per_division = insertion_probability,
    deletion_probability_per_target_per_division = deletion_probability,
    inter_target_dropout_probability = dropout_probability,
    inter_target_dropout_radius = dropout_radius,
    dropout_scope = "direct deletion cuts within one lineage branch segment",
    # How far a cut resects either side, as gamma shape/rate per side. NULL
    # keeps the legacy mean-one exponential, which will not reach the long
    # deletions a real array shows.
    deletion_extent = deletion_extent
  )

  for (cell_type in names(params$cell_type_dict$cell_type_params)) {
    cell_params <- params$cell_type_dict$cell_type_params[[cell_type]]
    cell_params$uninduced_editing_params$be_mutations_per_target_per_division <- 0
    cell_params$uninduced_editing_params$
      nuc_insertions_per_target_per_division <- 0
    cell_params$uninduced_editing_params$
      nuc_deletions_per_target_per_division <- 0
    cell_params$induced_editing_params$be_mutations_per_target_per_division <- 0
    cell_params$induced_editing_params$
      nuc_insertions_per_target_per_division <- insertion_probability
    cell_params$induced_editing_params$
      nuc_deletions_per_target_per_division <- deletion_probability
    params$cell_type_dict$cell_type_params[[cell_type]] <- cell_params
  }
  params
}

attach_flare_wt_crispr_model <- function(model, params) {
  configuration <- params$flare_wt_crispr_backend
  if (!isTRUE(configuration$enabled)) return(model)

  model$is_wt_crispr <- TRUE
  model$wt_crispr <- list(
    target_positions = sort(as.integer(names(model$nuc_target_classes))),
    insertion_probability_per_target_per_division =
      configuration$insertion_probability_per_target_per_division,
    deletion_probability_per_target_per_division =
      configuration$deletion_probability_per_target_per_division,
    inter_target_dropout_probability =
      configuration$inter_target_dropout_probability,
    inter_target_dropout_radius =
      configuration$inter_target_dropout_radius,
    dropout_scope = configuration$dropout_scope,
    deletion_extent = configuration$deletion_extent,
    barcode_length = length(model$barcode_sequence)
  )
  model
}

flare_wt_crispr_empty_events <- function() {
  data.frame(
    integration = integer(),
    position = integer(),
    event = character(),
    reference = character(),
    alternate = character(),
    allele = numeric(),
    event_time = numeric(),
    stringsAsFactors = FALSE
  )
}

mutate_flare_wt_crispr_segment <- function(profile,
                                            duration,
                                            rate_set,
                                            model,
                                            segment_start = 0) {
  if (duration <= 0 || length(rate_set$active_positions) == 0L) {
    return(list(profile = profile, events = flare_wt_crispr_empty_events()))
  }

  target_positions <- model$wt_crispr$target_positions
  dropout_probability <- model$wt_crispr$inter_target_dropout_probability
  dropout_radius <- model$wt_crispr$inter_target_dropout_radius
  event_rows <- list()
  event_count <- 0L

  append_event <- function(integration, position, event, alternate, allele,
                           event_time) {
    event_count <<- event_count + 1L
    event_rows[[event_count]] <<- data.frame(
      integration = as.integer(integration),
      position = as.integer(position),
      event = event,
      reference = model$barcode_sequence[position],
      alternate = alternate,
      allele = as.numeric(allele),
      event_time = as.numeric(event_time),
      stringsAsFactors = FALSE
    )
  }

  for (integration in seq_len(model$num_integrations)) {
    candidate_positions <- rate_set$active_positions[vapply(
      rate_set$active_positions,
      function(position) {
        physicell_barcode_value(profile, integration, position) == 0
      },
      logical(1)
    )]
    if (!length(candidate_positions)) next

    hazards <- rate_set$total_hazards[candidate_positions]
    mutation_draws <- stats::runif(length(candidate_positions))
    mutated <- mutation_draws < -expm1(-hazards * duration)
    if (!any(mutated)) next

    scheduled_positions <- candidate_positions[mutated]
    scheduled_times <- segment_start -
      log1p(-mutation_draws[mutated]) / hazards[mutated]
    event_order <- order(scheduled_times, scheduled_positions)
    scheduled_positions <- scheduled_positions[event_order]
    scheduled_times <- scheduled_times[event_order]
    direct_deletion_positions <- integer()

    for (event_index in seq_along(scheduled_positions)) {
      position <- scheduled_positions[event_index]
      event_time <- scheduled_times[event_index]

      # A prior inter-target deletion can remove a target whose cut was
      # independently scheduled later in this branch.
      if (physicell_barcode_value(profile, integration, position) != 0) next

      selected_event <- select_physicell_position_events(position, rate_set)
      decoded <- decode_physicell_barcode_events(selected_event)
      if (decoded$event == "deletion") {
        # A cut removes a contiguous run, not just its own base. Both extents
        # are drawn, so breakpoints scatter either side of the cut the way
        # resection puts them; the legacy port extended leftward only.
        extent <- draw_physicell_deletion_extent(
          model$wt_crispr$deletion_extent
        )
        applied <- apply_physicell_deletion(
          profile, integration, position,
          left_bases = extent[["left"]], right_bases = extent[["right"]],
          barcode_length = model$wt_crispr$barcode_length
        )
        profile <- applied$profile
        for (emptied in applied$positions) {
          append_event(
            integration, emptied,
            if (emptied == position) "deletion" else "resected_deletion",
            "-", -1, event_time
          )
        }
      } else {
        profile <- set_physicell_barcode_value(
          profile, integration, position, decoded$allele
        )
        append_event(
          integration,
          position,
          decoded$event,
          decoded$alternate,
          decoded$allele,
          event_time
        )
      }

      if (decoded$event != "deletion") next

      eligible_previous <- direct_deletion_positions[
        abs(direct_deletion_positions - position) <= dropout_radius
      ]
      for (previous_position in eligible_previous) {
        if (stats::runif(1L) >= dropout_probability) next

        bounds <- sort(c(previous_position, position))
        dropped_targets <- target_positions[
          target_positions >= bounds[[1L]] &
            target_positions <= bounds[[2L]]
        ]
        for (dropped_position in dropped_targets) {
          if (physicell_barcode_value(
              profile, integration, dropped_position
            ) == -1) {
            next
          }
          profile <- set_physicell_barcode_value(
            profile, integration, dropped_position, -1
          )
          append_event(
            integration,
            dropped_position,
            "inter_target_deletion",
            "-",
            -1,
            event_time
          )
        }
      }
      direct_deletion_positions <- c(direct_deletion_positions, position)
    }
  }

  events <- if (!event_count) {
    flare_wt_crispr_empty_events()
  } else {
    events <- do.call(rbind, event_rows)
    events[order(
      events$event_time,
      events$integration,
      events$position,
      events$event
    ), , drop = FALSE]
  }
  rownames(events) <- NULL
  list(profile = profile, events = events)
}

flare_wt_crispr_target_layout <- function(model) {
  target_positions <- model$wt_crispr$target_positions
  induced <- model$rate_sets$induced_editing_params$probabilities
  rate_value <- function(position, event) {
    value <- induced[[position]][event]
    if (!length(value) || is.na(value)) 0 else unname(value)
  }
  insertion_rates <- vapply(
    target_positions, rate_value, numeric(1), event = "insertion"
  )
  deletion_rates <- vapply(
    target_positions, rate_value, numeric(1), event = "deletion"
  )

  do.call(rbind, lapply(seq_len(model$num_integrations), function(integration) {
    data.frame(
      integration = integration,
      target_index = seq_along(target_positions),
      position = target_positions,
      edit_rate_class = unname(
        model$nuc_target_classes[as.character(target_positions)]
      ),
      reference = model$barcode_sequence[target_positions],
      alternate = "insertion|deletion",
      induced_insertion_probability_per_target_per_division = insertion_rates,
      induced_deletion_probability_per_target_per_division = deletion_rates,
      inter_target_dropout_probability =
        model$wt_crispr$inter_target_dropout_probability,
      inter_target_dropout_radius =
        model$wt_crispr$inter_target_dropout_radius,
      stringsAsFactors = FALSE
    )
  }))
}

flare_wt_crispr_logical_target_matrix <- function(state_matrix,
                                                   model,
                                                   integrations) {
  target_positions <- model$wt_crispr$target_positions
  barcode_length <- model$wt_crispr$barcode_length
  if (is.null(barcode_length)) {
    barcode_length <- max(target_positions)
  }

  # Each distinct deletion gets its own character state, keyed on the breakpoints
  # of the run that removed the target.
  #
  # Reading the target position alone -- which is what this did -- returns -1 for
  # every deletion, so a 156bp cut from 90 and a 183bp cut from 64 become the
  # same state and two cells that share neither ancestry nor event look
  # identical. The breakpoints are already in the profile, in the spacer
  # positions the deletion crossed; they only needed to be carried through. On
  # the HL60 comparison this is the difference between 28 and 311 distinguishable
  # whole-molecule outcomes, against 411 in the assay.
  #
  # States are assigned per target column from a shared dictionary, so the same
  # event in two cells maps to the same state and a shared derived character
  # still means shared ancestry.
  cell_count <- nrow(state_matrix)
  result <- matrix(0, cell_count, integrations * length(target_positions))
  dictionaries <- vector("list", length(target_positions))
  dictionaries[] <- list(new.env(hash = TRUE, parent = emptyenv()))
  next_state <- rep(0L, length(target_positions))

  for (integration in seq_len(integrations)) {
    wanted <- paste0("int_", integration, "_pos_", seq_len(barcode_length))
    columns <- match(wanted, colnames(state_matrix))
    if (all(is.na(columns))) {
      stop("WT-CRISPR logical targets are absent from the state matrix.")
    }
    present <- !is.na(columns)
    block <- as.matrix(state_matrix[, columns[present], drop = FALSE])
    coordinates <- seq_len(barcode_length)[present]
    target_index <- match(target_positions, coordinates)
    if (anyNA(target_index)) {
      stop("WT-CRISPR logical targets are absent from the state matrix.")
    }
    offset <- (integration - 1L) * length(target_positions)

    for (cell in seq_len(cell_count)) {
      values <- block[cell, ]
      states <- numeric(length(target_positions))
      deleted <- values == -1
      if (any(deleted)) {
        runs <- rle(deleted)
        ends <- cumsum(runs$lengths)
        starts <- ends - runs$lengths + 1L
        for (run in which(runs$values)) {
          covered <- which(target_index >= starts[run] & target_index <= ends[run])
          if (!length(covered)) next
          key <- sprintf("%dD+%d", runs$lengths[run], coordinates[starts[run]])
          for (target in covered) {
            store <- dictionaries[[target]]
            if (!exists(key, envir = store, inherits = FALSE)) {
              next_state[target] <- next_state[target] + 1L
              assign(key, next_state[target], envir = store)
            }
            states[target] <- get(key, envir = store, inherits = FALSE)
          }
        }
      }
      # An insertion is its own outcome and keeps a state of its own, offset past
      # the deletion states so the two cannot collide.
      inserted <- which(values[target_index] != 0 & values[target_index] != -1)
      for (target in inserted) {
        key <- sprintf("I%s", format(values[target_index][target]))
        store <- dictionaries[[target]]
        if (!exists(key, envir = store, inherits = FALSE)) {
          next_state[target] <- next_state[target] + 1L
          assign(key, next_state[target], envir = store)
        }
        states[target] <- get(key, envir = store, inherits = FALSE)
      }
      result[cell, offset + seq_along(target_positions)] <- states
    }
  }

  colnames(result) <- unlist(lapply(
    seq_len(integrations),
    function(integration) {
      paste0("int_", integration, "_target_", seq_along(target_positions))
    }
  ))
  rownames(result) <- rownames(state_matrix)
  # Downstream reaches into the @x slot, so this has to stay a sparse Matrix as
  # the subsetting version did; a base matrix fails there rather than here.
  Matrix::Matrix(result, sparse = TRUE)
}

annotate_flare_wt_crispr_manifest <- function(target_layout) {
  target_layout$logical_positions <- as.character(target_layout$position)
  target_layout$logical_state_encoding <- paste(
    "0=wild_type;-1=deletion;positive fractional allele=insertion;",
    "inter-target dropout is encoded as deletion"
  )
  target_layout
}

# Registration. Sourcing this file is now sufficient to dispatch the recorder:
# the simulator consults the registry at each of its four hook points, so there
# is nothing for a caller to install and nothing it can forget to install. The
# previous arrangement needed four wrappers assigned over simulator functions in
# another file, and configuring the parameters without them produced a run that
# looked right and silently used the native nuclease path throughout.
register_origin_recorder(
  "flare_wt_crispr",
  applies = function(model) isTRUE(model$is_wt_crispr),
  prepare = function(model, params) attach_flare_wt_crispr_model(model, params),
  mutate = mutate_flare_wt_crispr_segment,
  layout = flare_wt_crispr_target_layout,
  character_matrix = function(raw_alleles, binary_scores, model, integrations) {
    # FLARE_BINARY_CHARACTERS=1 restores the pre-allele encoding so the two can
    # be compared on identical simulated data.
    if (nzchar(Sys.getenv("FLARE_BINARY_CHARACTERS"))) {
      return(physicell_prime_character_matrix(raw_alleles, binary_scores, model))
    }
    flare_wt_crispr_logical_target_matrix(raw_alleles, model, integrations)
  }
)
