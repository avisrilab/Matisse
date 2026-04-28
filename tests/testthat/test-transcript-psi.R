# ---------------------------------------------------------------------------
# Tests for transcript mode construction and its internal helpers
# ---------------------------------------------------------------------------

# ---- .parse_ioe_files -------------------------------------------------------

test_that(".parse_ioe_files: returns correct event table structure", {
  f      <- make_ioe_file()
  events <- Matisse:::.parse_ioe_files(f)
  expect_s3_class(events, "data.frame")
  expect_equal(nrow(events), 2L)
  expect_true(all(c("event_id", "gene_id", "chr", "strand",
                    "event_type", "inclusion_features",
                    "exclusion_features") %in% colnames(events)))
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
  expect_equal(events$exclusion_features[1], "tx3")
  exc2 <- sort(strsplit(events$exclusion_features[2], ";")[[1]])
  expect_equal(exc2, c("tx6", "tx7"))
})

test_that(".parse_ioe_files: combines multiple IOE files with distinct events", {
  f1 <- make_ioe_file()
  # Second file with different event coordinates so no duplicates
  f2 <- tempfile(fileext = ".ioe")
  writeLines(c(
    "seqname\tgene_id\tinclusion_features\ttotal_transcripts",
    paste(c("chr2", "ENSG00000002;SE:chr2:100-200:300-400:+",
            "tx8,tx9", "tx8,tx9,tx10"), collapse = "\t")
  ), f2)
  events <- Matisse:::.parse_ioe_files(c(f1, f2))
  expect_equal(nrow(events), 3L)   # 2 from f1 + 1 from f2
})

test_that(".parse_ioe_files: deduplicates events shared across IOE files", {
  # Identical files -> all event IDs are duplicates; only unique ones kept
  f1 <- make_ioe_file()
  f2 <- make_ioe_file()
  expect_warning(
    Matisse:::.parse_ioe_files(c(f1, f2)),
    regexp = "duplicate"
  )
  events <- suppressWarnings(Matisse:::.parse_ioe_files(c(f1, f2)))
  expect_equal(nrow(events), 2L)   # 2 unique events, duplicates dropped
})

test_that(".parse_ioe_files: errors on malformed gene_id column", {
  bad <- tempfile(fileext = ".ioe")
  writeLines(c(
    "seqname\tgene_id\tinclusion_features\ttotal_transcripts",
    "chr1\tNO_SEMICOLON_HERE\ttx1\ttx1,tx2"
  ), bad)
  expect_error(Matisse:::.parse_ioe_files(bad), regexp = "malformed")
})

test_that(".parse_ioe_files: handles 5-column IOE format (SUPPA2 v2+)", {
  f <- tempfile(fileext = ".ioe")
  writeLines(c(
    "seqname\tgene_id\tevent_id\tinclusion_features\ttotal_transcripts",
    "chr1\tENSG00000001\tENSG00000001;SE:chr1:100-200:300-400:+\ttx1,tx2\ttx1,tx2,tx3"
  ), f)
  events <- Matisse:::.parse_ioe_files(f)
  expect_equal(nrow(events), 1L)
  expect_equal(events$gene_id,   "ENSG00000001")
  expect_equal(events$event_id,  "SE:chr1:100-200:300-400:+")
  expect_equal(events$inclusion_features, "tx1;tx2")
  expect_equal(events$exclusion_features, "tx3")
})

test_that(".parse_ioe_files: error message shows row number and bad value", {
  bad <- tempfile(fileext = ".ioe")
  writeLines(c(
    "seqname\tgene_id\tinclusion_features\ttotal_transcripts",
    "chr1\tENSG00000001;SE:chr1:100-200:300-400:+\ttx1\ttx1,tx2",
    "chr2\tNO_SEMICOLON_HERE\ttx3\ttx3,tx4"
  ), bad)
  err <- tryCatch(Matisse:::.parse_ioe_files(bad), error = conditionMessage)
  expect_match(err, "row 2")
  expect_match(err, "NO_SEMICOLON_HERE")
  expect_match(err, basename(bad))
})

# ---- .aggregate_transcript_counts ------------------------------------------
# Post-Phase 6: .parse_ioe_files emits inclusion_features / exclusion_features
# directly, so .aggregate_transcript_counts can consume the parser output
# without renaming. The previous .rename_ioe_for_aggregate bridge is gone.

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

# ---- CreateMatisseObject (transcript mode) ----------------------------------

test_that("CreateMatisseObject (transcript mode): returns a MatisseObject", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_from_transcripts()
  expect_s4_class(obj, "MatisseObject")
})

test_that("CreateMatisseObject (transcript mode): input.mode is 'transcript'", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_from_transcripts()
  expect_equal(obj@input.mode, "transcript")
})

test_that("CreateMatisseObject (transcript mode): 'isoform' assay stored in Seurat", {
  skip_if_not_installed("Seurat")
  obj      <- make_matisse_from_transcripts()
  iso_assay <- GetSeurat(obj)[["isoform"]]
  expect_false(is.null(iso_assay))
  expect_true(inherits(iso_assay, "Assay5"))
})

test_that("CreateMatisseObject (transcript mode): nCount_isoform written at construction", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  obj    <- CreateMatisseObject(
    seurat            = seu,
    transcript_counts = tx_mat,
    events         = make_ioe_file(),
    verbose           = FALSE
  )
  expect_true("nCount_isoform"   %in% colnames(MatisseMeta(obj)))
  expect_true("nFeature_isoform" %in% colnames(MatisseMeta(obj)))
})

test_that("CreateMatisseObject (transcript mode): isoform assay holds transcripts x cells", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_from_transcripts()
  tx  <- SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]], "counts")
  expect_equal(ncol(tx), 10L)   # cells
  expect_equal(nrow(tx), 8L)    # transcripts
})

test_that("CalculatePSI (transcript mode): 'psi' Assay5 stored in Seurat", {
  skip_if_not_installed("Seurat")
  obj       <- make_matisse_from_transcripts()
  psi_assay <- GetSeurat(obj)[["psi"]]
  expect_false(is.null(psi_assay))
  expect_true(inherits(psi_assay, "Assay5"))
})

test_that("CalculatePSI (transcript mode): GetPSI returns cells x events", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_from_transcripts()
  psi <- GetPSI(obj)
  expect_equal(nrow(psi), 10L)
  expect_equal(ncol(psi), 2L)
})

test_that("CalculatePSI (transcript mode): PSI values are in [0,1] or NA", {
  skip_if_not_installed("Seurat")
  obj    <- make_matisse_from_transcripts()
  vals   <- as.numeric(GetPSI(obj))
  finite <- vals[!is.na(vals)]
  expect_true(all(finite >= 0 & finite <= 1))
})

test_that("CalculatePSI (transcript mode): nPercent_isoform written to metadata", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_from_transcripts()
  expect_true("nPercent_isoform" %in% colnames(MatisseMeta(obj)))
  vals <- MatisseMeta(obj)$nPercent_isoform
  expect_true(all(vals >= 0 & vals <= 100))
})

test_that("CalculatePSI (transcript mode): event annotation populated in PSI assay meta.features", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_from_transcripts()
  mf  <- GetSeurat(obj)[["psi"]][[]]
  expect_equal(nrow(mf), 2L)
  expect_true(all(c("gene_id", "chr", "strand",
                    "event_type", "inclusion_features",
                    "exclusion_features") %in% colnames(mf)))
})

test_that("CalculatePSI (transcript mode): inclusion + exclusion sums to total for covered entries", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()
  obj    <- CreateMatisseObject(
    seurat = seu, transcript_counts = tx_mat,
    events = f, verbose = FALSE)
  obj    <- CalculatePSI(obj, min_coverage = 0L, verbose = FALSE)

  inc     <- as.matrix(Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "counts")))
  exc     <- as.matrix(Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "exclusion")))
  psi_mat <- as.matrix(GetPSI(obj))

  covered <- !is.na(psi_mat) & (inc + exc) > 0
  expect_true(all(abs(inc[covered] / (inc[covered] + exc[covered]) -
                        psi_mat[covered]) < 1e-9))
})

test_that("CreateMatisseObject (transcript mode): errors if seurat is wrong type", {
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()
  expect_error(
    CreateMatisseObject(
      seurat = list(), transcript_counts = tx_mat, events = f),
    regexp = "must be a Seurat object")
})

test_that("CreateMatisseObject (transcript mode): errors if no cell overlap", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  f      <- make_ioe_file()
  bad_mat <- Matrix::Matrix(
    matrix(1L, nrow = 3, ncol = 5,
           dimnames = list(paste0("tx", 1:3), paste0("X", 1:5))),
    sparse = TRUE)
  expect_error(
    CreateMatisseObject(
      seurat = seu, transcript_counts = bad_mat, events = f),
    regexp = "No cell barcodes overlap")
})

test_that("CreateMatisseObject (transcript mode): errors if IOE file missing", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  expect_error(
    CreateMatisseObject(
      seurat = seu, transcript_counts = tx_mat,
      events = "/nonexistent/path.ioe"),
    regexp = "not found")
})

test_that("CreateMatisseObject (transcript mode): warns on partial cell overlap", {
  skip_if_not_installed("Seurat")
  seu        <- make_seurat()
  tx_partial <- make_transcript_counts()[, paste0("Cell", 1:9), drop = FALSE]
  f          <- make_ioe_file()
  expect_warning(
    CreateMatisseObject(
      seurat = seu, transcript_counts = tx_partial,
      events = f, verbose = FALSE),
    regexp = "9/10")
})

test_that("CreateMatisseObject (transcript mode): passes validObject check", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  f      <- make_ioe_file()
  obj    <- CreateMatisseObject(
    seurat = seu, transcript_counts = tx_mat,
    events = f, verbose = FALSE)
  expect_no_error(methods::validObject(obj))
})
