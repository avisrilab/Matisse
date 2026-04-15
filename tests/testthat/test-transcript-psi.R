# ---------------------------------------------------------------------------
# Tests for CreateMatisseObjectFromTranscripts and its internal helpers
# ---------------------------------------------------------------------------

# ---- .parse_ioe_files -------------------------------------------------------

test_that(".parse_ioe_files: returns correct event table structure", {
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)

  expect_s3_class(events, "data.frame")
  expect_equal(nrow(events), 2L)
  expect_true(all(c("event_id", "gene_id", "chr", "strand",
                    "event_type", "inclusion_transcripts",
                    "exclusion_transcripts") %in% colnames(events)))
})

test_that(".parse_ioe_files: parses event_id and gene_id correctly", {
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)

  expect_equal(events$gene_id[1],    "ENSG00000001")
  expect_equal(events$event_id[1],   "SE:chr1:100-200:300-400:+")
  expect_equal(events$event_type[1], "SE")
  expect_equal(events$strand[1],     "+")
  expect_equal(events$chr[1],        "chr1")
})

test_that(".parse_ioe_files: exclusion = total minus inclusion", {
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)

  # Event 1: total = tx1,tx2,tx3; inclusion = tx1,tx2 → exclusion = tx3
  expect_equal(events$exclusion_transcripts[1], "tx3")

  # Event 2: total = tx4,tx5,tx6,tx7; inclusion = tx4,tx5 → exclusion = tx6;tx7
  exc2 <- sort(strsplit(events$exclusion_transcripts[2], ";")[[1]])
  expect_equal(exc2, c("tx6", "tx7"))
})

test_that(".parse_ioe_files: combines multiple IOE files", {
  f1 <- make_ioe_file()
  f2 <- make_ioe_file()   # identical content → 4 rows total
  events <- Matisse:::.parse_ioe_files(c(f1, f2))
  expect_equal(nrow(events), 4L)
})

test_that(".parse_ioe_files: errors on malformed gene_id column", {
  bad <- tempfile(fileext = ".ioe")
  writeLines(c(
    "seqname\tgene_id\tinclusion_transcripts\ttotal_transcripts",
    "chr1\tNO_SEMICOLON_HERE\ttx1\ttx1,tx2"
  ), bad)
  expect_error(Matisse:::.parse_ioe_files(bad), regexp = "malformed")
})

# ---- .aggregate_transcript_counts ------------------------------------------

test_that(".aggregate_transcript_counts: output dimensions are cells x events", {
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)
  tx_mat <- make_transcript_counts()
  cells  <- colnames(tx_mat)

  res <- Matisse:::.aggregate_transcript_counts(
    tx_counts = tx_mat, events = events,
    min_coverage = 1L, cells = cells)

  expect_equal(dim(res$psi),       c(10L, 2L))
  expect_equal(dim(res$inclusion), c(10L, 2L))
  expect_equal(dim(res$exclusion), c(10L, 2L))
})

test_that(".aggregate_transcript_counts: PSI = 1 when only inclusion reads", {
  cells  <- "Cell1"
  # tx1=10, tx2=5, tx3=0 → PSI = 15/15 = 1
  mat    <- Matrix::Matrix(matrix(c(10, 5, 0, 0, 0, 0, 0, 0),
                                   nrow = 8, ncol = 1,
                                   dimnames = list(paste0("tx", 1:8), cells)),
                            sparse = TRUE)
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)

  res <- Matisse:::.aggregate_transcript_counts(
    mat, events, min_coverage = 1L, cells = cells)

  expect_equal(as.numeric(res$psi["Cell1", "SE:chr1:100-200:300-400:+"]), 1.0)
})

test_that(".aggregate_transcript_counts: PSI = 0 when only exclusion reads", {
  cells  <- "Cell1"
  # tx1=0, tx2=0, tx3=10 → PSI = 0/10 = 0
  mat    <- Matrix::Matrix(matrix(c(0, 0, 10, 0, 0, 0, 0, 0),
                                   nrow = 8, ncol = 1,
                                   dimnames = list(paste0("tx", 1:8), cells)),
                            sparse = TRUE)
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)

  res <- Matisse:::.aggregate_transcript_counts(
    mat, events, min_coverage = 1L, cells = cells)

  expect_equal(as.numeric(res$psi["Cell1", "SE:chr1:100-200:300-400:+"]), 0.0)
})

test_that(".aggregate_transcript_counts: PSI = 0.5 with equal counts", {
  cells <- "Cell1"
  # tx1=5, tx2=5 (inc=10), tx3=10 (exc=10) → PSI = 10/20 = 0.5
  mat   <- Matrix::Matrix(matrix(c(5, 5, 10, 0, 0, 0, 0, 0),
                                  nrow = 8, ncol = 1,
                                  dimnames = list(paste0("tx", 1:8), cells)),
                           sparse = TRUE)
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)

  res <- Matisse:::.aggregate_transcript_counts(
    mat, events, min_coverage = 1L, cells = cells)

  expect_equal(as.numeric(res$psi["Cell1", "SE:chr1:100-200:300-400:+"]), 0.5)
})

test_that(".aggregate_transcript_counts: low-coverage entries become NA", {
  cells <- "Cell1"
  # tx1=1, tx2=1, tx3=1 → total = 3 < min_coverage 5 → NA
  mat   <- Matrix::Matrix(matrix(c(1, 1, 1, 0, 0, 0, 0, 0),
                                  nrow = 8, ncol = 1,
                                  dimnames = list(paste0("tx", 1:8), cells)),
                           sparse = TRUE)
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)

  res <- Matisse:::.aggregate_transcript_counts(
    mat, events, min_coverage = 5L, cells = cells)

  expect_true(is.na(res$psi["Cell1", "SE:chr1:100-200:300-400:+"]))
})

test_that(".aggregate_transcript_counts: missing transcripts treated as zero", {
  cells <- "Cell1"
  # Only tx3 present (exclusion); inclusion transcripts absent → PSI = 0
  mat   <- Matrix::Matrix(matrix(c(0, 0, 8, 0, 0, 0, 0, 0),
                                  nrow = 8, ncol = 1,
                                  dimnames = list(paste0("tx", 1:8), cells)),
                           sparse = TRUE)
  # Remove inclusion transcripts from the matrix entirely
  mat_partial <- mat[c("tx3"), , drop = FALSE]

  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)[1, ]  # only event 1

  res <- Matisse:::.aggregate_transcript_counts(
    mat_partial, events, min_coverage = 1L, cells = cells)

  expect_equal(as.numeric(res$psi[1, 1]), 0.0)
})

# ---- CreateMatisseObjectFromTranscripts -------------------------------------

test_that("CreateMatisseObjectFromTranscripts: returns a MatisseObject", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()

  obj <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 1L, verbose = FALSE)

  expect_s4_class(obj, "MatisseObject")
})

test_that("CreateMatisseObjectFromTranscripts: PSI dimensions are cells x events", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()

  obj <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 1L, verbose = FALSE)

  expect_equal(nrow(obj@psi), 10L)
  expect_equal(ncol(obj@psi), 2L)
})

test_that("CreateMatisseObjectFromTranscripts: PSI values are in [0,1] or NA", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()

  obj    <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 1L, verbose = FALSE)

  vals   <- as.numeric(obj@psi)
  finite <- vals[!is.na(vals)]
  expect_true(all(finite >= 0 & finite <= 1))
})

test_that("CreateMatisseObjectFromTranscripts: event_data is populated", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()

  obj <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 1L, verbose = FALSE)

  ed <- GetEventData(obj)
  expect_equal(nrow(ed), 2L)
  expect_true(all(c("event_id", "gene_id", "chr", "strand",
                    "event_type", "inclusion_junctions",
                    "exclusion_junctions") %in% colnames(ed)))
})

test_that("CreateMatisseObjectFromTranscripts: junction_counts slot is NULL", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()

  obj <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 1L, verbose = FALSE)

  expect_null(obj@junction_counts)
})

test_that("CreateMatisseObjectFromTranscripts: inclusion + exclusion = total", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()

  obj    <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 0L, verbose = FALSE)

  inc <- as.matrix(obj@inclusion_counts)
  exc <- as.matrix(obj@exclusion_counts)
  psi <- as.matrix(obj@psi)

  # Where PSI is non-NA, inc/(inc+exc) should equal PSI
  covered <- !is.na(psi) & (inc + exc) > 0
  expect_true(all(abs(inc[covered] / (inc[covered] + exc[covered]) -
                        psi[covered]) < 1e-9))
})

test_that("CreateMatisseObjectFromTranscripts: errors if seurat is wrong type", {
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()
  expect_error(
    CreateMatisseObjectFromTranscripts(
      seurat = list(), transcript_counts = tx_mat, ioe_files = f),
    regexp = "must be a Seurat object")
})

test_that("CreateMatisseObjectFromTranscripts: errors if no cell overlap", {
  seu    <- make_seurat()
  f      <- make_ioe_file()
  # Matrix with completely different cell names
  bad_mat <- Matrix::Matrix(
    matrix(1L, nrow = 3, ncol = 5,
           dimnames = list(paste0("tx", 1:3), paste0("X", 1:5))),
    sparse = TRUE)
  expect_error(
    CreateMatisseObjectFromTranscripts(
      seurat = seu, transcript_counts = bad_mat, ioe_files = f),
    regexp = "No cell barcodes overlap")
})

test_that("CreateMatisseObjectFromTranscripts: errors if IOE file missing", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  expect_error(
    CreateMatisseObjectFromTranscripts(
      seurat = seu, transcript_counts = tx_mat,
      ioe_files = "/nonexistent/path.ioe"),
    regexp = "not found")
})

test_that("CreateMatisseObjectFromTranscripts: warns on partial cell overlap", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()   # 10 cells: Cell1–Cell10
  tx_mat <- make_transcript_counts()

  # Drop Cell10 from transcript matrix
  tx_partial <- tx_mat[, colnames(tx_mat) != "Cell10", drop = FALSE]
  f <- make_ioe_file()

  expect_warning(
    CreateMatisseObjectFromTranscripts(
      seurat = seu, transcript_counts = tx_partial,
      ioe_files = f, min_coverage = 1L, verbose = FALSE),
    regexp = "9/10")
})

test_that("CreateMatisseObjectFromTranscripts: passes validObject check", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()

  obj <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 1L, verbose = FALSE)

  expect_no_error(methods::validObject(obj))
})

test_that("CreateMatisseObjectFromTranscripts: downstream QC works", {
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()

  obj <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 1L, verbose = FALSE)

  obj <- ComputeIsoformQC(obj, verbose = FALSE)
  expect_true("mean_psi" %in% colnames(MatisseMeta(obj)))
})
