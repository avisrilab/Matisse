# ---------------------------------------------------------------------------
# Tests for QC functions (ComputeIsoformQC, FilterCells, FilterEvents)
# ---------------------------------------------------------------------------

test_that("ComputeIsoformQC: adds expected columns to isoform_metadata", {
  obj       <- make_matisse_object()
  obj       <- CalculatePSI(obj, verbose = FALSE)
  obj       <- ComputeIsoformQC(obj, verbose = FALSE)
  meta_cols <- colnames(MatisseMeta(obj))
  expect_true("n_junctions_detected"  %in% meta_cols)
  expect_true("total_junction_reads"  %in% meta_cols)
  expect_true("n_events_covered"      %in% meta_cols)
  expect_true("pct_events_covered"    %in% meta_cols)
  expect_true("mean_psi"              %in% meta_cols)
})

test_that("ComputeIsoformQC: n_junctions_detected is non-negative integer", {
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, verbose = FALSE)
  obj  <- ComputeIsoformQC(obj, verbose = FALSE)
  vals <- MatisseMeta(obj)$n_junctions_detected
  expect_true(all(vals >= 0L))
  expect_type(vals, "integer")
})

test_that("ComputeIsoformQC: total_junction_reads >= n_junctions_detected", {
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, verbose = FALSE)
  obj  <- ComputeIsoformQC(obj, verbose = FALSE)
  meta <- MatisseMeta(obj)
  expect_true(all(meta$total_junction_reads >= meta$n_junctions_detected))
})

test_that("ComputeIsoformQC: pct_events_covered is in [0, 100]", {
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, verbose = FALSE)
  obj  <- ComputeIsoformQC(obj, verbose = FALSE)
  vals <- MatisseMeta(obj)$pct_events_covered
  expect_true(all(vals >= 0 & vals <= 100))
})

test_that("ComputeIsoformQC: mean_psi is in [0, 1] or NA", {
  obj    <- make_matisse_object()
  obj    <- CalculatePSI(obj, verbose = FALSE)
  obj    <- ComputeIsoformQC(obj, verbose = FALSE)
  vals   <- MatisseMeta(obj)$mean_psi
  finite <- vals[!is.nan(vals) & !is.na(vals)]
  expect_true(all(finite >= 0 & finite <= 1))
})

test_that("ComputeIsoformQC: warns when PSI assay is absent", {
  obj <- make_matisse_object()  # PSI not yet calculated
  expect_warning(
    ComputeIsoformQC(obj, verbose = TRUE),
    regexp = "PSI matrix is NULL"
  )
})

test_that("ComputeIsoformQC: row count matches number of cells", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  obj <- ComputeIsoformQC(obj, verbose = FALSE)
  expect_equal(nrow(MatisseMeta(obj)), .n_cells(obj))
})

test_that("ComputeIsoformQC: works on transcript-based object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  obj <- ComputeIsoformQC(obj, verbose = FALSE)
  expect_true("mean_psi" %in% colnames(MatisseMeta(obj)))
  # No junction metrics expected since junction_counts is NULL
  expect_false("n_junctions_detected" %in% colnames(MatisseMeta(obj)))
})

# ---------------------------------------------------------------------------
# FilterCells
# ---------------------------------------------------------------------------

test_that("FilterCells: removes cells below min_junctions threshold", {
  obj    <- make_matisse_object()
  obj    <- CalculatePSI(obj, verbose = FALSE)
  obj    <- ComputeIsoformQC(obj, verbose = FALSE)
  thresh <- max(MatisseMeta(obj)$n_junctions_detected)
  sub    <- FilterCells(obj, min_junctions = thresh, verbose = FALSE)
  expect_true(.n_cells(sub) <= .n_cells(obj))
  if (.n_cells(sub) > 0) {
    expect_true(all(MatisseMeta(sub)$n_junctions_detected >= thresh))
  }
})

test_that("FilterCells: with permissive threshold keeps all cells", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  obj <- ComputeIsoformQC(obj, verbose = FALSE)
  sub <- FilterCells(obj, min_junctions = 0L, verbose = FALSE)
  expect_equal(.n_cells(sub), .n_cells(obj))
})

test_that("FilterCells: warns for unknown QC column", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  expect_warning(
    FilterCells(obj,
      custom_filters = list(n_junctions_detected = c(1, NA)),
      verbose = FALSE),
    regexp = "not found in isoform_metadata"
  )
})

test_that("FilterCells: custom_filters with upper bound works", {
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, verbose = FALSE)
  obj  <- ComputeIsoformQC(obj, verbose = FALSE)
  med  <- stats::median(MatisseMeta(obj)$total_junction_reads)
  sub  <- FilterCells(obj,
    custom_filters = list(total_junction_reads = c(NA, med)),
    verbose = FALSE)
  if (.n_cells(sub) > 0) {
    expect_true(all(MatisseMeta(sub)$total_junction_reads <= med))
  }
})

test_that("FilterCells: 'psi' assay cells are also subsetted", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  obj <- ComputeIsoformQC(obj, verbose = FALSE)
  sub <- FilterCells(obj, min_junctions = 0L, max_junctions = 3L,
                     verbose = FALSE)
  # After filtering, PSI assay cell count should match
  psi <- GetPSI(sub)
  expect_equal(nrow(psi), .n_cells(sub))
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

  # Manually set first event to constant PSI = 0.5 via SetPSI
  psi_cx       <- as.matrix(GetPSI(obj))
  psi_cx[, 1]  <- 0.5
  psi_sp       <- Matrix::Matrix(psi_cx, sparse = TRUE)
  obj          <- SetPSI(obj, psi_sp)

  sub <- FilterEvents(obj, min_cells_covered = 0L,
                      min_psi_variance = 0.001,
                      verbose = FALSE)
  expect_true(.n_events(sub) < .n_events(obj))
})
