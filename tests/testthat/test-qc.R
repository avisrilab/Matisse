# ---------------------------------------------------------------------------
# Tests for QC columns and filtering (FilterCells, FilterEvents)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# QC columns written at construction and by CalculatePSI
# ---------------------------------------------------------------------------

test_that("CreateMatisseObject writes nCount_isoform and nFeature_isoform", {
  obj       <- make_matisse_object()
  meta_cols <- colnames(MatisseMeta(obj))
  expect_true("nCount_isoform"   %in% meta_cols)
  expect_true("nFeature_isoform" %in% meta_cols)
})

test_that("nCount_isoform is a non-negative integer", {
  obj  <- make_matisse_object()
  vals <- MatisseMeta(obj)$nCount_isoform
  expect_type(vals, "integer")
  expect_true(all(vals >= 0L))
})

test_that("nFeature_isoform <= nCount_isoform per cell", {
  obj  <- make_matisse_object()
  meta <- MatisseMeta(obj)
  expect_true(all(meta$nFeature_isoform <= meta$nCount_isoform))
})

test_that("CalculatePSI writes nPercent_isoform", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, min_coverage = 1L, verbose = FALSE)
  expect_true("nPercent_isoform" %in% colnames(MatisseMeta(obj)))
})

test_that("nPercent_isoform is in [0, 100]", {
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, min_coverage = 1L, verbose = FALSE)
  vals <- MatisseMeta(obj)$nPercent_isoform
  expect_true(all(vals >= 0 & vals <= 100))
})

test_that("QC columns row count matches number of cells", {
  obj <- make_matisse_object()
  expect_equal(nrow(MatisseMeta(obj)), .n_cells(obj))
})

test_that("transcript mode: nCount_isoform and nFeature_isoform written at construction", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  obj <- CreateMatisseObject(
    seurat            = seu,
    transcript_counts = tx_mat,
    ioe_files         = make_ioe_file(),
    verbose           = FALSE
  )
  meta_cols <- colnames(MatisseMeta(obj))
  expect_true("nCount_isoform"   %in% meta_cols)
  expect_true("nFeature_isoform" %in% meta_cols)
})

test_that("transcript mode: CalculatePSI writes nPercent_isoform", {
  obj <- make_matisse_from_transcripts()  # already called CalculatePSI
  expect_true("nPercent_isoform" %in% colnames(MatisseMeta(obj)))
})

# ---------------------------------------------------------------------------
# FilterCells
# ---------------------------------------------------------------------------

test_that("FilterCells: removes cells below min_features_isoform threshold", {
  obj    <- make_matisse_object()
  thresh <- max(MatisseMeta(obj)$nFeature_isoform)
  sub    <- FilterCells(obj, min_features_isoform = thresh, verbose = FALSE)
  expect_true(.n_cells(sub) <= .n_cells(obj))
  if (.n_cells(sub) > 0L) {
    expect_true(all(MatisseMeta(sub)$nFeature_isoform >= thresh))
  }
})

test_that("FilterCells: permissive threshold keeps all cells", {
  obj <- make_matisse_object()
  sub <- FilterCells(obj, min_features_isoform = 0L, verbose = FALSE)
  expect_equal(.n_cells(sub), .n_cells(obj))
})

test_that("FilterCells: max_counts_isoform upper bound works", {
  obj <- make_matisse_object()
  med <- stats::median(MatisseMeta(obj)$nCount_isoform)
  sub <- FilterCells(obj, max_counts_isoform = med, verbose = FALSE)
  if (.n_cells(sub) > 0L) {
    expect_true(all(MatisseMeta(sub)$nCount_isoform <= med))
  }
})

test_that("FilterCells: min_pct_isoform filters on nPercent_isoform", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, min_coverage = 1L, verbose = FALSE)
  sub <- FilterCells(obj, min_pct_isoform = 0, verbose = FALSE)
  expect_equal(.n_cells(sub), .n_cells(obj))
})

test_that("FilterCells: warns for unknown QC column in custom_filters", {
  obj <- make_matisse_object()
  expect_warning(
    FilterCells(obj,
      custom_filters = list(nonexistent_col = c(1, NA)),
      verbose = FALSE),
    regexp = "not found in metadata"
  )
})

test_that("FilterCells: 'psi' assay cells are also subsetted", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  sub <- FilterCells(obj, min_features_isoform = 0L, max_features_isoform = 3L,
                     verbose = FALSE)
  psi <- GetPSI(sub)
  if (!is.null(psi)) expect_equal(nrow(psi), .n_cells(sub))
})

test_that("FilterCells: no warning in transcript mode", {
  obj <- make_matisse_from_transcripts()
  expect_no_warning(FilterCells(obj, verbose = FALSE))
})

# ---------------------------------------------------------------------------
# FilterEvents
# ---------------------------------------------------------------------------

test_that("FilterEvents: errors if 'psi' assay is absent", {
  obj <- make_matisse_object()
  expect_error(FilterEvents(obj), regexp = "PSI matrix is NULL")
})

test_that("FilterEvents: with min_cells_covered = 0 keeps all events", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  sub <- FilterEvents(obj, min_cells_covered = 0L, verbose = FALSE)
  expect_equal(.n_events(sub), .n_events(obj))
})

test_that("FilterEvents: with impossibly high threshold drops all events", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  sub <- FilterEvents(obj, min_cells_covered = .n_cells(obj) + 1L,
                      verbose = FALSE)
  expect_equal(.n_events(sub), 0L)
})

test_that("FilterEvents: variance filter removes low-variance events", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, min_coverage = 1L, verbose = FALSE)

  # Set first event to constant PSI = 0.5 via SetPSI
  psi_cx       <- as.matrix(GetPSI(obj))
  psi_cx[, 1]  <- 0.5
  psi_sp       <- Matrix::Matrix(psi_cx, sparse = TRUE)
  obj          <- SetPSI(obj, psi_sp)

  sub <- FilterEvents(obj, min_cells_covered = 0L,
                      min_psi_variance = 0.001,
                      verbose = FALSE)
  expect_true(.n_events(sub) < .n_events(obj))
})
