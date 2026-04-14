# ---------------------------------------------------------------------------
# Tests for QC functions
# ---------------------------------------------------------------------------

test_that("ComputeIsoformQC: adds expected columns to isoform_metadata", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  obj <- ComputeIsoformQC(obj, verbose = FALSE)
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
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, verbose = FALSE)
  obj  <- ComputeIsoformQC(obj, verbose = FALSE)
  vals <- MatisseMeta(obj)$mean_psi
  finite <- vals[!is.nan(vals) & !is.na(vals)]
  expect_true(all(finite >= 0 & finite <= 1))
})

test_that("ComputeIsoformQC: warns when PSI is NULL", {
  obj <- make_matisse_object()  # PSI not yet calculated
  expect_warning(
    ComputeIsoformQC(obj, verbose = TRUE),
    regexp = "PSI matrix is NULL"
  )
})

test_that("ComputeIsoformQC: row count matches number of cells", {
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, verbose = FALSE)
  obj  <- ComputeIsoformQC(obj, verbose = FALSE)
  expect_equal(nrow(MatisseMeta(obj)), .n_cells(obj))
})

# ---------------------------------------------------------------------------
# FilterCells
# ---------------------------------------------------------------------------

test_that("FilterCells: removes cells below min_junctions threshold", {
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, verbose = FALSE)
  obj  <- ComputeIsoformQC(obj, verbose = FALSE)
  meta <- MatisseMeta(obj)

  # Set threshold just above minimum to guarantee removal
  thresh <- max(meta$n_junctions_detected)
  sub    <- FilterCells(obj, min_junctions = thresh, verbose = FALSE)

  expect_true(.n_cells(sub) <= .n_cells(obj))
  if (.n_cells(sub) > 0) {
    expect_true(all(
      MatisseMeta(sub)$n_junctions_detected >= thresh))
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
  # isoform_metadata has no QC columns yet without ComputeIsoformQC
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
  meta <- MatisseMeta(obj)

  # Keep only cells with total_junction_reads <= median
  med   <- stats::median(meta$total_junction_reads)
  sub   <- FilterCells(obj,
    custom_filters = list(total_junction_reads = c(NA, med)),
    verbose = FALSE)
  if (.n_cells(sub) > 0) {
    expect_true(all(
      MatisseMeta(sub)$total_junction_reads <= med))
  }
})

# ---------------------------------------------------------------------------
# FilterEvents
# ---------------------------------------------------------------------------

test_that("FilterEvents: errors if PSI is NULL", {
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
  # n_cells + 1 as threshold: no event can be covered in more cells than exist
  sub <- FilterEvents(obj, min_cells_covered = .n_cells(obj) + 1L,
                      verbose = FALSE)
  expect_equal(.n_events(sub), 0L)
})

test_that("FilterEvents: variance filter removes low-variance events", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)

  # Manually set one event to constant PSI = 0.5 (zero variance)
  psi  <- as.matrix(obj@psi)
  psi[, 1] <- 0.5
  obj@psi  <- Matrix::Matrix(psi, sparse = TRUE)

  sub <- FilterEvents(obj, min_cells_covered = 0L,
                      min_psi_variance = 0.001,
                      verbose = FALSE)
  expect_true(.n_events(sub) < .n_events(obj))
})
