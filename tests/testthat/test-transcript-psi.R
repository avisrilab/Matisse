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
  expect_equal(events$exclusion_transcripts[1], "tx3")
  exc2 <- sort(strsplit(events$exclusion_transcripts[2], ";")[[1]])
  expect_equal(exc2, c("tx6", "tx7"))
})

test_that(".parse_ioe_files: combines multiple IOE files", {
  f1     <- make_ioe_file()
  f2     <- make_ioe_file()
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

test_that(".parse_ioe_files: handles 5-column IOE format (SUPPA2 v2+)", {
  f <- tempfile(fileext = ".ioe")
  writeLines(c(
    "seqname\tgene_id\tevent_id\tinclusion_transcripts\ttotal_transcripts",
    "chr1\tENSG00000001\tENSG00000001;SE:chr1:100-200:300-400:+\ttx1,tx2\ttx1,tx2,tx3"
  ), f)
  events <- Matisse:::.parse_ioe_files(f)
  expect_equal(nrow(events), 1L)
  expect_equal(events$gene_id,   "ENSG00000001")
  expect_equal(events$event_id,  "SE:chr1:100-200:300-400:+")
  expect_equal(events$inclusion_transcripts, "tx1;tx2")
  expect_equal(events$exclusion_transcripts, "tx3")
})

test_that(".parse_ioe_files: error message shows row number and bad value", {
  bad <- tempfile(fileext = ".ioe")
  writeLines(c(
    "seqname\tgene_id\tinclusion_transcripts\ttotal_transcripts",
    "chr1\tENSG00000001;SE:chr1:100-200:300-400:+\ttx1\ttx1,tx2",
    "chr2\tNO_SEMICOLON_HERE\ttx3\ttx3,tx4"
  ), bad)
  err <- tryCatch(Matisse:::.parse_ioe_files(bad), error = conditionMessage)
  expect_match(err, "row 2")
  expect_match(err, "NO_SEMICOLON_HERE")
  expect_match(err, basename(bad))
})

# ---- .aggregate_transcript_counts ------------------------------------------

test_that(".aggregate_transcript_counts: output dimensions are cells x events", {
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)
  tx_mat <- make_transcript_counts()
  cells  <- colnames(tx_mat)
  res <- Matisse:::.aggregate_transcript_counts(
    tx_counts = tx_mat, events = events, min_coverage = 1L, cells = cells)
  expect_equal(dim(res$psi),       c(10L, 2L))
  expect_equal(dim(res$inclusion), c(10L, 2L))
  expect_equal(dim(res$exclusion), c(10L, 2L))
})

test_that(".aggregate_transcript_counts: PSI = 1 when only inclusion reads", {
  cells <- "Cell1"
  mat   <- Matrix::Matrix(matrix(c(10, 5, 0, 0, 0, 0, 0, 0),
                                  nrow = 8, ncol = 1,
                                  dimnames = list(paste0("tx", 1:8), cells)),
                           sparse = TRUE)
  events <- Matisse:::.parse_ioe_files(make_ioe_file())
  res    <- Matisse:::.aggregate_transcript_counts(
    mat, events, min_coverage = 1L, cells = cells)
  expect_equal(as.numeric(res$psi["Cell1", "SE:chr1:100-200:300-400:+"]), 1.0)
})

test_that(".aggregate_transcript_counts: PSI = 0 when only exclusion reads", {
  cells <- "Cell1"
  mat   <- Matrix::Matrix(matrix(c(0, 0, 10, 0, 0, 0, 0, 0),
                                  nrow = 8, ncol = 1,
                                  dimnames = list(paste0("tx", 1:8), cells)),
                           sparse = TRUE)
  events <- Matisse:::.parse_ioe_files(make_ioe_file())
  res    <- Matisse:::.aggregate_transcript_counts(
    mat, events, min_coverage = 1L, cells = cells)
  expect_equal(as.numeric(res$psi["Cell1", "SE:chr1:100-200:300-400:+"]), 0.0)
})

test_that(".aggregate_transcript_counts: PSI = 0.5 with equal counts", {
  cells <- "Cell1"
  mat   <- Matrix::Matrix(matrix(c(5, 5, 10, 0, 0, 0, 0, 0),
                                  nrow = 8, ncol = 1,
                                  dimnames = list(paste0("tx", 1:8), cells)),
                           sparse = TRUE)
  events <- Matisse:::.parse_ioe_files(make_ioe_file())
  res    <- Matisse:::.aggregate_transcript_counts(
    mat, events, min_coverage = 1L, cells = cells)
  expect_equal(as.numeric(res$psi["Cell1", "SE:chr1:100-200:300-400:+"]), 0.5)
})

test_that(".aggregate_transcript_counts: low-coverage entries become NA", {
  cells <- "Cell1"
  mat   <- Matrix::Matrix(matrix(c(1, 1, 1, 0, 0, 0, 0, 0),
                                  nrow = 8, ncol = 1,
                                  dimnames = list(paste0("tx", 1:8), cells)),
                           sparse = TRUE)
  events <- Matisse:::.parse_ioe_files(make_ioe_file())
  res    <- Matisse:::.aggregate_transcript_counts(
    mat, events, min_coverage = 5L, cells = cells)
  expect_true(is.na(res$psi["Cell1", "SE:chr1:100-200:300-400:+"]))
})

test_that(".aggregate_transcript_counts: missing transcripts treated as zero", {
  cells       <- "Cell1"
  mat_partial <- Matrix::Matrix(matrix(c(8), nrow = 1, ncol = 1,
                                        dimnames = list("tx3", cells)),
                                 sparse = TRUE)
  events <- Matisse:::.parse_ioe_files(make_ioe_file())[1, ]
  res    <- Matisse:::.aggregate_transcript_counts(
    mat_partial, events, min_coverage = 1L, cells = cells)
  expect_equal(as.numeric(res$psi[1, 1]), 0.0)
})

# ---- CreateMatisseObjectFromTranscripts -------------------------------------

test_that("CreateMatisseObjectFromTranscripts: returns a MatisseObject", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  expect_s4_class(obj, "MatisseObject")
})

test_that("CreateMatisseObjectFromTranscripts: 'transcript' assay stored in Seurat", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj      <- make_matisse_from_transcripts()
  tx_assay <- GetSeurat(obj)[["transcript"]]
  expect_false(is.null(tx_assay))
})

test_that("CreateMatisseObjectFromTranscripts: GetTranscriptCounts returns transcripts x cells", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  tx  <- GetTranscriptCounts(obj)
  expect_equal(ncol(tx), 10L)   # cells
  expect_equal(nrow(tx), 8L)    # transcripts
})

test_that("CreateMatisseObjectFromTranscripts: 'psi' Assay5 stored in Seurat", {
  skip_if_not_installed("Seurat")
  obj       <- make_matisse_from_transcripts()
  psi_assay <- GetSeurat(obj)[["psi"]]
  expect_false(is.null(psi_assay))
  expect_true(inherits(psi_assay, "Assay5"))
})

test_that("CreateMatisseObjectFromTranscripts: GetPSI returns cells x events", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  psi <- GetPSI(obj)
  expect_equal(nrow(psi), 10L)
  expect_equal(ncol(psi), 2L)
})

test_that("CreateMatisseObjectFromTranscripts: PSI values are in [0,1] or NA", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj    <- make_matisse_from_transcripts()
  vals   <- as.numeric(GetPSI(obj))
  finite <- vals[!is.na(vals)]
  expect_true(all(finite >= 0 & finite <= 1))
})

test_that("CreateMatisseObjectFromTranscripts: event_data is populated", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  ed  <- GetEventData(obj)
  expect_equal(nrow(ed), 2L)
  expect_true(all(c("event_id", "gene_id", "chr", "strand",
                    "event_type", "inclusion_junctions",
                    "exclusion_junctions") %in% colnames(ed)))
})

test_that("CreateMatisseObjectFromTranscripts: junction_counts slot is NULL", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  expect_null(GetJunctionCounts(obj))
})

test_that("CreateMatisseObjectFromTranscripts: inclusion + exclusion sums to total for covered entries", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()
  obj    <- CreateMatisseObjectFromTranscripts(
    seurat = seu, transcript_counts = tx_mat,
    ioe_files = f, min_coverage = 0L, verbose = FALSE)

  inc     <- as.matrix(GetInclusionCounts(obj))
  exc     <- as.matrix(GetExclusionCounts(obj))
  psi_mat <- as.matrix(GetPSI(obj))

  covered <- !is.na(psi_mat) & (inc + exc) > 0
  expect_true(all(abs(inc[covered] / (inc[covered] + exc[covered]) -
                        psi_mat[covered]) < 1e-9))
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
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  f      <- make_ioe_file()
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
  skip_if_not_installed("Seurat")
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
  skip_if_not_installed("Signac")
  seu        <- make_seurat()
  tx_partial <- make_transcript_counts()[, paste0("Cell", 1:9), drop = FALSE]
  f          <- make_ioe_file()
  expect_warning(
    CreateMatisseObjectFromTranscripts(
      seurat = seu, transcript_counts = tx_partial,
      ioe_files = f, min_coverage = 1L, verbose = FALSE),
    regexp = "9/10")
})

test_that("CreateMatisseObjectFromTranscripts: passes validObject check", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  expect_no_error(methods::validObject(obj))
})

test_that("CreateMatisseObjectFromTranscripts: downstream QC works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  obj <- ComputeIsoformQC(obj, verbose = FALSE)
  expect_true("mean_psi" %in% colnames(MatisseMeta(obj)))
})
