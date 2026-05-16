#' @include MatisseObject-class.R
#' @include constructors.R
NULL

# ---------------------------------------------------------------------------
# SUPPA2 events -> junction-ID event table (short-read / junction mode)
# ---------------------------------------------------------------------------

#' Build a junction-ID event table from SUPPA2 events
#'
#' SUPPA2 \code{.ioe} files describe each splicing event with \emph{transcript}
#' IDs, but junction mode needs the inclusion/exclusion features to be
#' \emph{junction} IDs that match a STARsolo SJ matrix. The genomic
#' coordinates of every defining junction are, however, fully encoded in the
#' SUPPA2 \code{event_id} (e.g.
#' \code{SE:chr1:804222-804776:804966-807217:+}). This function parses that
#' grammar per event type and emits the
#' \code{chr-start-end-strand} junction IDs Matisse expects, returning the
#' \code{data.frame} that \code{\link{CreateMatisseObject}} accepts via
#' \code{events}. This is a deterministic coordinate adapter, not a heuristic
#' event caller.
#'
#' \strong{Coordinate convention.} SUPPA2 reports exon-boundary positions
#' (last base of the donor exon, first base of the acceptor exon). STARsolo SJ
#' reports the first/last \emph{intronic} base. The two differ by a fixed
#' offset (intron = \code{[donor + 1, acceptor - 1]}). When
#' \code{junction_universe} is supplied and \code{intron_offset = "auto"}, the
#' offset is calibrated by maximising overlap with the observed STARsolo
#' junctions; otherwise the standard \code{c(1L, 1L)} is used.
#'
#' \strong{Supported event types.} \code{SE}, \code{A3}, \code{A5},
#' \code{MX}, \code{AF}, \code{AL}. \code{RI} (retained intron) is \emph{not}
#' junction-quantifiable - intron retention has no supporting split junction -
#' and is dropped with a warning if requested.
#'
#' @param ioe_files Character vector of SUPPA2 \code{.ioe} file path(s).
#' @param event_types Character vector of event types to keep. Default:
#'   \code{c("SE", "A3", "A5", "MX", "AF", "AL")}.
#' @param junction_universe Optional character vector of observed junction IDs
#'   (e.g. \code{colnames(ReadSTARsoloSJ(...))}). Used to auto-calibrate
#'   \code{intron_offset} and to report match rate. Default: \code{NULL}.
#' @param intron_offset Either \code{"auto"} (calibrate against
#'   \code{junction_universe}; falls back to \code{c(1L, 1L)} if absent) or an
#'   integer vector \code{c(left, right)} added to the donor / subtracted from
#'   the acceptor exon-boundary coordinate. Default: \code{"auto"}.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return A \code{data.frame} with columns \code{event_id}, \code{gene_id},
#'   \code{chr}, \code{strand}, \code{event_type}, \code{inclusion_features},
#'   \code{exclusion_features} (the last two semicolon-separated junction IDs).
#'
#' @seealso \code{\link{ReadSTARsoloSJ}}, \code{\link{CreateMatisseObject}}.
#'
#' @examples
#' \dontrun{
#' jxn    <- ReadSTARsoloSJ("Solo.out/SJ/raw", cells = colnames(seu))
#' events <- BuildJunctionEvents(
#'   list.files("events", "\\.ioe$", full.names = TRUE),
#'   junction_universe = colnames(jxn)
#' )
#' obj <- CreateMatisseObject(seu, junction_counts = jxn, events = events)
#' }
#'
#' @export
BuildJunctionEvents <- function(ioe_files,
                                event_types = c("SE", "A3", "A5",
                                                "MX", "AF", "AL"),
                                junction_universe = NULL,
                                intron_offset = "auto",
                                verbose = TRUE) {
  if (!is.character(ioe_files) || length(ioe_files) == 0L) {
    rlang::abort("`ioe_files` must be a non-empty character vector of paths.")
  }
  missing_files <- ioe_files[!file.exists(ioe_files)]
  if (length(missing_files) > 0L) {
    rlang::abort(paste0("IOE file(s) not found: ",
                        paste(missing_files, collapse = ", ")))
  }

  if (verbose) cli::cli_alert_info("Parsing {length(ioe_files)} IOE file(s)...")
  ev <- .parse_ioe_files(ioe_files)         # reuses tested SUPPA2 parser

  unsupported <- setdiff(event_types, .matisse_junction_event_types())
  if (length(unsupported) > 0L) {
    rlang::warn(paste0(
      "Ignoring event type(s) not junction-quantifiable: ",
      paste(unsupported, collapse = ", "),
      ". Supported: ",
      paste(.matisse_junction_event_types(), collapse = ", "), "."))
    event_types <- intersect(event_types, .matisse_junction_event_types())
  }
  keep <- ev$event_type %in% event_types
  if ("RI" %in% ev$event_type[!keep] && verbose) {
    cli::cli_alert_warning(paste0(
      "Dropped RI events: intron retention has no supporting split ",
      "junction, so it cannot be scored from STARsolo SJ counts."))
  }
  ev <- ev[keep, , drop = FALSE]
  if (nrow(ev) == 0L) {
    rlang::abort("No events of a supported type remained after filtering.")
  }

  # Parse each event_id into inclusion / exclusion donor-acceptor coordinate
  # pairs (still in SUPPA2 exon-boundary space).
  pairs <- .suppa_event_pairs(ev$event_id, ev$event_type)

  off <- .resolve_intron_offset(intron_offset, pairs, ev$chr, ev$strand,
                                junction_universe, verbose)

  ev$inclusion_features <- .pairs_to_junction_ids(
    pairs$inc, ev$chr, ev$strand, off)
  ev$exclusion_features <- .pairs_to_junction_ids(
    pairs$exc, ev$chr, ev$strand, off)

  if (!is.null(junction_universe) && verbose) {
    rate <- .junction_match_rate(ev, junction_universe)
    cli::cli_alert_info(paste0(
      "Junction match rate against STARsolo universe: ",
      sprintf("%.1f%%", 100 * rate),
      " of defining junctions observed."))
  }

  out <- data.frame(
    event_id           = ev$event_id,
    gene_id            = ev$gene_id,
    chr                = ev$chr,
    strand             = ev$strand,
    event_type         = ev$event_type,
    inclusion_features = ev$inclusion_features,
    exclusion_features = ev$exclusion_features,
    stringsAsFactors   = FALSE,
    row.names          = NULL
  )
  if (verbose) {
    cli::cli_alert_success(paste0(
      "Built {nrow(out)} junction-mode events ",
      "({paste(sort(unique(out$event_type)), collapse = '/')})."))
  }
  out
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.matisse_junction_event_types <- function() {
  c("SE", "A3", "A5", "MX", "AF", "AL")
}

# Parse SUPPA2 event_ids into per-event inclusion / exclusion coordinate
# pairs. Returns list(inc = list(<event> -> list(c(donor, acceptor), ...)),
# exc = same). Coordinates are SUPPA2 exon-boundary positions.
.suppa_event_pairs <- function(event_ids, event_types) {
  n   <- length(event_ids)
  inc <- vector("list", n)
  exc <- vector("list", n)
  for (i in seq_len(n)) {
    toks   <- strsplit(event_ids[i], ":", fixed = TRUE)[[1L]]
    # toks: [type, chr, <coord blocks...>, strand]
    blocks <- toks[3:(length(toks) - 1L)]
    je     <- .suppa_blocks_to_junctions(event_types[i], blocks)
    inc[[i]] <- je$inc
    exc[[i]] <- je$exc
  }
  list(inc = inc, exc = exc)
}

# Map one event's coordinate blocks to inclusion / exclusion junctions.
# Each junction is an integer vector c(donor_boundary, acceptor_boundary)
# in SUPPA2 exon-boundary space (5'->3' on the genome, not strand-oriented;
# SUPPA2 already writes coordinates left-to-right on the + reference).
.suppa_blocks_to_junctions <- function(type, blocks) {
  is_pair <- grepl("-", blocks, fixed = TRUE)
  prs <- lapply(blocks[is_pair], function(b) {
    as.integer(strsplit(b, "-", fixed = TRUE)[[1L]])
  })
  if (type == "SE") {
    # SE:e1-s2:e2-s3  -> inc {e1-s2, e2-s3}, exc {e1-s3}
    e1 <- prs[[1L]][1L]; s2 <- prs[[1L]][2L]
    e2 <- prs[[2L]][1L]; s3 <- prs[[2L]][2L]
    list(inc = list(c(e1, s2), c(e2, s3)), exc = list(c(e1, s3)))
  } else if (type == "MX") {
    # MX:e1-s2:e2-s3:e1-s4:e4-s3 -> inc {pair1,pair2}, exc {pair3,pair4}
    list(inc = list(prs[[1L]], prs[[2L]]),
         exc = list(prs[[3L]], prs[[4L]]))
  } else {
    # A3 / A5 / AF / AL: each isoform contributes one discriminating
    # junction; first pair is the inclusion (SUPPA2 "alternative") form,
    # second the exclusion form.
    list(inc = list(prs[[1L]]), exc = list(prs[[2L]]))
  }
}

# Resolve / calibrate the exon-boundary -> intron offset.
.resolve_intron_offset <- function(intron_offset, pairs, chr, strand,
                                    junction_universe, verbose) {
  if (is.numeric(intron_offset)) {
    if (length(intron_offset) != 2L) {
      rlang::abort("`intron_offset` must be length-2 c(left, right).")
    }
    return(as.integer(intron_offset))
  }
  if (!identical(intron_offset, "auto")) {
    rlang::abort('`intron_offset` must be "auto" or an integer c(left, right).')
  }
  if (is.null(junction_universe)) return(c(1L, 1L))

  candidates <- list(c(1L, 1L), c(0L, 0L), c(1L, 0L),
                     c(0L, 1L), c(-1L, -1L))
  uni  <- unique(junction_universe)
  idx  <- seq_len(min(length(chr), 3000L))
  flat <- lapply(idx, function(i) c(pairs$inc[[i]], pairs$exc[[i]]))
  rates <- vapply(candidates, function(cand) {
    ids <- unlist(lapply(idx, function(k) {
      vapply(flat[[k]], function(p)
        .junction_id(chr[k], p[1L] + cand[1L],
                     p[2L] - cand[2L], strand[k]),
        character(1))
    }), use.names = FALSE)
    mean(ids %in% uni)
  }, numeric(1))
  best      <- candidates[[which.max(rates)]]
  best_rate <- max(rates)
  if (verbose) {
    spread <- paste(vapply(seq_along(candidates), function(i)
      sprintf("c(%d,%d)=%.1f%%", candidates[[i]][1L], candidates[[i]][2L],
              100 * rates[i]), character(1)), collapse = "  ")
    cli::cli_alert_info(paste0("Offset calibration spread: ", spread))
    # A correct convention should win decisively over the alternatives;
    # a near-tie means the annotation/genome build likely mismatches.
    runner_up <- max(rates[-which.max(rates)])
    if (best_rate > 0 && runner_up / best_rate > 0.75) {
      cli::cli_alert_warning(paste0(
        "Best offset is not decisive (runner-up within 75%); check that the ",
        "SUPPA2 annotation matches the STARsolo genome build."))
    }
    cli::cli_alert_info(paste0(
      "Calibrated intron offset c(", best[1L], ", ", best[2L], ") ",
      "(", sprintf("%.1f%%", 100 * best_rate),
      " of sampled junctions matched)."))
  }
  best
}

.junction_id <- function(chr, start, end, strand) {
  paste(chr, start, end, strand, sep = "-")
}

# Convert a list-of-events of coordinate pairs to ";"-joined junction IDs.
.pairs_to_junction_ids <- function(pair_lists, chr, strand, off) {
  vapply(seq_along(pair_lists), function(i) {
    ps <- pair_lists[[i]]
    if (length(ps) == 0L) return("")
    ids <- vapply(ps, function(p)
      .junction_id(chr[i], p[1L] + off[1L], p[2L] - off[2L], strand[i]),
      character(1))
    paste(unique(ids), collapse = ";")
  }, character(1))
}

.junction_match_rate <- function(ev, junction_universe) {
  uni <- unique(junction_universe)
  ids <- unlist(strsplit(
    c(ev$inclusion_features, ev$exclusion_features), ";", fixed = TRUE),
    use.names = FALSE)
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0L) return(0)
  mean(ids %in% uni)
}
