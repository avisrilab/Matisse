# ---------------------------------------------------------------------------
# Tests for MergeMatisse
# ---------------------------------------------------------------------------

# Helper: two independent MatisseObjects with PSI calculated
make_pair <- function() {
  obj1 <- CalculatePSI(make_matisse_object(), min_coverage = 1L, verbose = FALSE)
  obj2 <- CalculatePSI(make_matisse_object(), min_coverage = 1L, verbose = FALSE)
  list(obj1 = obj1, obj2 = obj2)
}

# ---- Basic structure -------------------------------------------------------

test_that("MergeMatisse: merged cell count equals the sum of both objects", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(.n_cells(merged), .n_cells(p$obj1) + .n_cells(p$obj2))
})

test_that("MergeMatisse: PSI row count equals sum of both objects", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(nrow(GetPSI(merged)), .n_cells(p$obj1) + .n_cells(p$obj2))
})

test_that("MergeMatisse: PSI event (column) count is preserved", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(ncol(GetPSI(merged)), .n_events(p$obj1))
})

.jxn_counts <- function(obj) {
  Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]], "counts"))
}

test_that("MergeMatisse: junction_counts row count equals sum of both objects", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(
    nrow(.jxn_counts(merged)),
    nrow(.jxn_counts(p$obj1)) + nrow(.jxn_counts(p$obj2))
  )
})

test_that("MergeMatisse: junction_counts column count (junctions) is preserved", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(ncol(.jxn_counts(merged)), ncol(.jxn_counts(p$obj1)))
})

# ---- Cell naming -----------------------------------------------------------

test_that("MergeMatisse: default add_cell_ids = c('x','y') prefixes cell names", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  cells  <- .get_cells(merged)
  expect_true(any(startsWith(cells, "x_")))
  expect_true(any(startsWith(cells, "y_")))
})

test_that("MergeMatisse: custom add_cell_ids are applied", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, add_cell_ids = c("s1", "s2"), verbose = FALSE)
  cells  <- .get_cells(merged)
  expect_true(any(startsWith(cells, "s1_")))
  expect_true(any(startsWith(cells, "s2_")))
})

# ---- Data integrity --------------------------------------------------------

test_that("MergeMatisse: PSI values remain in [0, 1] or NA", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  psi    <- GetPSI(merged)
  vals   <- as.numeric(psi)
  finite <- vals[!is.na(vals)]
  expect_true(all(finite >= 0 & finite <= 1))
})

test_that("MergeMatisse: event annotation is carried forward from x", {
  # Post-P1, event annotation lives in the PSI assay's meta.features.
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(GetSeurat(merged)[["psi"]][[]],
               GetSeurat(p$obj1)[["psi"]][[]])
})

test_that("MergeMatisse: per-junction coordinates carry through (P7)", {
  # After P7, per-junction coords live in the isoform assay's meta.features
  # (auto-parsed from junction IDs). They should ride through the Seurat
  # merge automatically.
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  iso_mf <- GetSeurat(merged)[["isoform"]][[]]
  expect_equal(rownames(iso_mf),
               rownames(GetSeurat(p$obj1)[["isoform"]][[]]))
  expect_true(all(c("chr", "start", "end", "strand") %in% colnames(iso_mf)))
  expect_true(all(!is.na(iso_mf$chr)))
})

test_that("MergeMatisse: returns a valid MatisseObject", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_s4_class(merged, "MatisseObject")
  expect_no_error(methods::validObject(merged))
})

# ---- Mismatched events -----------------------------------------------------
# Coverage gap: testing partial event-set overlap on merge requires shared
# >= 2 events (Assay5 minimum-features constraint) plus 3+ events overall.
# The current fixtures have only 2 events, so the previous tests used a
# single-event subset which is no longer representable. Re-cover the
# warn-on-mismatch + intersect-keep behaviour after extending the fixtures.

# ---- Error handling --------------------------------------------------------

test_that("MergeMatisse: errors when input.mode values differ (P10)", {
  skip_if_not_installed("Seurat")
  obj_jxn <- make_matisse_object()
  obj_tx  <- make_matisse_from_transcripts()
  expect_error(
    MergeMatisse(obj_jxn, obj_tx, verbose = FALSE),
    regexp = "different input.mode"
  )
})

test_that("MergeMatisse: errors when add_cell_ids is wrong length", {
  p <- make_pair()
  expect_error(MergeMatisse(p$obj1, p$obj2, add_cell_ids = "x"))
  expect_error(MergeMatisse(p$obj1, p$obj2, add_cell_ids = c("a", "b", "c")))
})

test_that("MergeMatisse: errors if x is not a MatisseObject", {
  p <- make_pair()
  expect_error(MergeMatisse("string", p$obj2), regexp = "MatisseObjects")
})

test_that("MergeMatisse: errors if y is not a MatisseObject", {
  p <- make_pair()
  expect_error(MergeMatisse(p$obj1, list()), regexp = "MatisseObjects")
})
