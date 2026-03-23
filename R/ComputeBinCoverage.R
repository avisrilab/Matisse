#' Compute per-cell bin support from a BAM file.
#'
#' This function computes the unnormalized per-cell bin-level inputs used by the
#' downstream transcript-resolution and plotting-normalization workflow.
#'
#' For each cell x bin combination, it returns:
#' 1. A signal value defined as the sum over reads of (overlap width / bin width)
#' 2. A read-support value defined as the number of distinct reads contributing
#'    support to the bin
#'
#' These outputs are intentionally left unnormalized here. Group scaling,
#' global scaling, and smoothing are applied later in the plotting workflow.
#'
#' @param bam Character(1). Path to a BAM file.
#' @param bins A `GenomicRanges::GRanges` object of genomic bins. Must contain a
#'   `bin_id` metadata column.
#' @param region Optional `GenomicRanges::GRanges` object specifying the genomic
#'   region to read from the BAM. If `NULL`, the range of `bins` is used.
#' @param barcode.field Integer(1). Position of the cell barcode after splitting
#'   each read name by `qname.sep`.
#' @param qname.sep Character(1). Delimiter used to split read names.
#'
#' @return A list containing:
#'   \describe{
#'     \item{signal_matrix}{Sparse cell-by-bin matrix of summed fractional bin support.}
#'     \item{read_matrix}{Sparse cell-by-bin matrix of distinct read counts contributing to each bin.}
#'     \item{cells}{Character vector of cell barcodes.}
#'     \item{bins}{Input `GRanges` of bins.}
#'   }
#' @export
ComputeBinCoverage <- function(
    bam,
    bins,
    region = NULL,
    barcode.field = 2,
    qname.sep = "_"
) {

  # ---- Helper checks ----
  .check_path <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1) file path.", call. = FALSE)
    }
    if (!file.exists(x)) {
      stop("File not found: ", x, call. = FALSE)
    }
  }

  .check_granges <- function(x, nm) {
    if (!methods::is(x, "GRanges")) {
      stop(nm, " must be a GenomicRanges::GRanges object.", call. = FALSE)
    }
  }

  .check_positive_integer <- function(x, nm) {
    if (length(x) != 1L || is.na(x) || x < 1L || x != as.integer(x)) {
      stop(nm, " must be a positive integer(1).", call. = FALSE)
    }
  }

  .check_scalar_character <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1).", call. = FALSE)
    }
  }

  # ---- Validate inputs ----
  .check_path(bam, "bam")
  .check_granges(bins, "bins")
  .check_positive_integer(barcode.field, "barcode.field")
  .check_scalar_character(qname.sep, "qname.sep")

  if (!is.null(region)) {
    .check_granges(region, "region")
    if (length(region) < 1L) {
      stop("region must contain at least one range when supplied.", call. = FALSE)
    }
  }

  if (length(bins) < 1L) {
    stop("bins must contain at least one range.", call. = FALSE)
  }

  bins.mcols <- S4Vectors::mcols(bins)
  if (!("bin_id" %in% colnames(bins.mcols))) {
    stop("bins must contain a 'bin_id' metadata column.", call. = FALSE)
  }

  if (anyNA(bins.mcols$bin_id) || any(!nzchar(as.character(bins.mcols$bin_id)))) {
    stop("bins$bin_id contains missing or empty values.", call. = FALSE)
  }

  # ---- use bins span if region not supplied ----
  if (is.null(region)) {
    region <- range(bins)
  }

  # ---- create BAM reading parameters ----
  param <- Rsamtools::ScanBamParam(
    which = region,
    what = c("qname")
  )

  # ---- read alignments from BAM ----
  gal <- GenomicAlignments::readGAlignments(
    file = bam,
    use.names = TRUE,
    param = param
  )

  # ---- extract read names and parse cell barcodes ----
  qnames <- names(gal)

  if (length(qnames) > 0L) {
    if (anyNA(qnames) || any(!nzchar(qnames))) {
      stop("Read names in BAM contain missing or empty values.", call. = FALSE)
    }

    qname.parts <- strsplit(qnames, split = qname.sep, fixed = TRUE)

    bad.qnames <- vapply(
      qname.parts,
      function(x) length(x) < barcode.field,
      logical(1)
    )

    if (any(bad.qnames)) {
      stop(
        "At least one read name does not contain enough fields for barcode.field = ",
        barcode.field, " when split by qname.sep = '", qname.sep, "'.",
        call. = FALSE
      )
    }
  } else {
    qname.parts <- list()
  }

  cell.barcodes <- vapply(
    qname.parts,
    function(x) x[[barcode.field]],
    character(1)
  )

  if (length(cell.barcodes) > 0L) {
    if (anyNA(cell.barcodes) || any(!nzchar(cell.barcodes))) {
      stop("Parsed cell barcodes contain missing or empty values.", call. = FALSE)
    }
  }

  all.cells <- sort(unique(cell.barcodes))
  all.bins <- as.character(S4Vectors::mcols(bins)$bin_id)

  # ---- early return if no reads ----
  if (length(gal) == 0L) {
    zero.signal <- Matrix::Matrix(
      0,
      nrow = length(all.cells),
      ncol = length(bins),
      sparse = TRUE,
      dimnames = list(all.cells, all.bins)
    )

    zero.reads <- Matrix::Matrix(
      0,
      nrow = length(all.cells),
      ncol = length(bins),
      sparse = TRUE,
      dimnames = list(all.cells, all.bins)
    )

    return(list(
      signal_matrix = zero.signal,
      read_matrix = zero.reads,
      cells = all.cells,
      bins = bins
    ))
  }

  # ---- convert alignments into genomic blocks ----
  gal.grl <- GenomicAlignments::grglist(gal, order.as.in.query = TRUE)

  # ---- flatten blocks into a GRanges object ----
  block.gr <- unlist(gal.grl, use.names = FALSE)
  read.id <- rep(seq_along(gal.grl), lengths(gal.grl))
  block.cell <- cell.barcodes[read.id]

  S4Vectors::mcols(block.gr)$read_id <- read.id
  S4Vectors::mcols(block.gr)$cell <- block.cell

  # ---- find which blocks overlap which bins ----
  hits <- GenomicRanges::findOverlaps(block.gr, bins, ignore.strand = FALSE)

  if (length(hits) == 0L) {
    zero.signal <- Matrix::Matrix(
      0,
      nrow = length(all.cells),
      ncol = length(bins),
      sparse = TRUE,
      dimnames = list(all.cells, all.bins)
    )

    zero.reads <- Matrix::Matrix(
      0,
      nrow = length(all.cells),
      ncol = length(bins),
      sparse = TRUE,
      dimnames = list(all.cells, all.bins)
    )

    return(list(
      signal_matrix = zero.signal,
      read_matrix = zero.reads,
      cells = all.cells,
      bins = bins
    ))
  }

  # ---- compute overlap widths and fractional bin contributions ----
  block.hit.gr <- block.gr[S4Vectors::queryHits(hits)]
  bin.hit.gr <- bins[S4Vectors::subjectHits(hits)]

  overlap.width <- IRanges::width(IRanges::pintersect(block.hit.gr, bin.hit.gr))
  bin.width <- BiocGenerics::width(bin.hit.gr)

  overlap.df <- data.frame(
    cell = S4Vectors::mcols(block.hit.gr)$cell,
    read_id = S4Vectors::mcols(block.hit.gr)$read_id,
    bin_id = as.character(S4Vectors::mcols(bin.hit.gr)$bin_id),
    overlap_width = overlap.width,
    bin_width = bin.width,
    stringsAsFactors = FALSE
  )

  overlap.df$fractional_value <- overlap.df$overlap_width / overlap.df$bin_width

  # ---- collapse to distinct read x cell x bin first ----
  # This avoids double-counting the same read within a bin if multiple blocks hit it.
  read.bin.df <- stats::aggregate(
    fractional_value ~ cell + read_id + bin_id,
    data = overlap.df,
    FUN = sum
  )

  # ---- cap per-read contribution at 1 ----
  # A read covering an entire bin contributes 1.
  read.bin.df$fractional_value <- pmin(read.bin.df$fractional_value, 1)

  # ---- aggregate signal per cell x bin ----
  signal.sum <- stats::aggregate(
    fractional_value ~ cell + bin_id,
    data = read.bin.df,
    FUN = sum
  )

  # ---- aggregate distinct read support per cell x bin ----
  read.sum <- stats::aggregate(
    read_id ~ cell + bin_id,
    data = read.bin.df,
    FUN = length
  )
  colnames(read.sum)[colnames(read.sum) == "read_id"] <- "read_count"

  # ---- build sparse signal matrix ----
  i.signal <- match(signal.sum$cell, all.cells)
  j.signal <- match(signal.sum$bin_id, all.bins)
  x.signal <- signal.sum$fractional_value

  signal.mat <- Matrix::sparseMatrix(
    i = i.signal,
    j = j.signal,
    x = x.signal,
    dims = c(length(all.cells), length(all.bins)),
    dimnames = list(all.cells, all.bins)
  )

  # ---- build sparse read-support matrix ----
  i.read <- match(read.sum$cell, all.cells)
  j.read <- match(read.sum$bin_id, all.bins)
  x.read <- read.sum$read_count

  read.mat <- Matrix::sparseMatrix(
    i = i.read,
    j = j.read,
    x = x.read,
    dims = c(length(all.cells), length(all.bins)),
    dimnames = list(all.cells, all.bins)
  )

  # ---- return ----
  list(
    signal_matrix = signal.mat,
    read_matrix = read.mat,
    cells = all.cells,
    bins = bins
  )
}
