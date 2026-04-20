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

test_that("MergeMatisse: junction_counts row count equals sum of both objects", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(
    nrow(GetJunctionCounts(merged)),
    nrow(GetJunctionCounts(p$obj1)) + nrow(GetJunctionCounts(p$obj2))
  )
})

test_that("MergeMatisse: junction_counts column count (junctions) is preserved", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(ncol(GetJunctionCounts(merged)), ncol(GetJunctionCounts(p$obj1)))
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

test_that("MergeMatisse: event_data is carried forward from x", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(GetEventData(merged), GetEventData(p$obj1))
})

test_that("MergeMatisse: junction_data is carried forward from x", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_equal(GetJunctionData(merged), GetJunctionData(p$obj1))
})

test_that("MergeMatisse: returns a valid MatisseObject", {
  p      <- make_pair()
  merged <- MergeMatisse(p$obj1, p$obj2, verbose = FALSE)
  expect_s4_class(merged, "MatisseObject")
  expect_no_error(methods::validObject(merged))
})

# ---- Mismatched events -----------------------------------------------------

test_that("MergeMatisse: warns when PSI column sets differ", {
  p    <- make_pair()
  # Subset obj2 to only one event
  obj2_sub <- p$obj2[, "SE_gene1_e2"]
  expect_warning(
    MergeMatisse(p$obj1, obj2_sub, verbose = FALSE),
    regexp = "shared events"
  )
})

test_that("MergeMatisse: result contains only shared events after mismatch", {
  p        <- make_pair()
  obj2_sub <- p$obj2[, "SE_gene1_e2"]
  merged   <- suppressWarnings(
    MergeMatisse(p$obj1, obj2_sub, verbose = FALSE)
  )
  expect_equal(.n_events(merged), 1L)
  expect_equal(colnames(GetPSI(merged)), "SE_gene1_e2")
})

# ---- Error handling --------------------------------------------------------

test_that("MergeMatisse: errors if x is not a MatisseObject", {
  p <- make_pair()
  expect_error(MergeMatisse("string", p$obj2), regexp = "MatisseObjects")
})

test_that("MergeMatisse: errors if y is not a MatisseObject", {
  p <- make_pair()
  expect_error(MergeMatisse(p$obj1, list()), regexp = "MatisseObjects")
})
