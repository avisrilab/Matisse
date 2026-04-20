# ---------------------------------------------------------------------------
# Tests for PSI calculation (CalculatePSI, SummarizePSI, internal helpers)
# ---------------------------------------------------------------------------

# ---- .parse_junction_list --------------------------------------------------

test_that(".parse_junction_list: parses semicolon-separated string", {
  result <- Matisse:::.parse_junction_list("jxn1;jxn2;jxn3")
  expect_equal(result, c("jxn1", "jxn2", "jxn3"))
})

test_that(".parse_junction_list: handles single junction", {
  result <- Matisse:::.parse_junction_list("jxn1")
  expect_equal(result, "jxn1")
})

test_that(".parse_junction_list: returns empty vector for NA input", {
  expect_equal(Matisse:::.parse_junction_list(NA_character_), character(0))
})

test_that(".parse_junction_list: returns empty vector for empty string", {
  expect_equal(Matisse:::.parse_junction_list(""), character(0))
})

test_that(".parse_junction_list: trims whitespace", {
  result <- Matisse:::.parse_junction_list(" jxn1 ; jxn2 ")
  expect_equal(result, c("jxn1", "jxn2"))
})

# ---- PSI on a raw matrix ---------------------------------------------------

test_that("CalculatePSI (matrix): returns a matrix of correct dimensions", {
  jxn_mat <- make_junction_counts()
  events  <- make_event_data()
  result  <- CalculatePSI(jxn_mat, events, min_coverage = 1L, verbose = FALSE)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(10L, 2L))
})

test_that("CalculatePSI (matrix): values are in [0, 1] or NA", {
  jxn_mat     <- make_junction_counts()
  events      <- make_event_data()
  result      <- CalculatePSI(jxn_mat, events, min_coverage = 1L, verbose = FALSE)
  finite_vals <- result[!is.na(result)]
  expect_true(all(finite_vals >= 0))
  expect_true(all(finite_vals <= 1))
})

test_that("CalculatePSI (matrix): low-coverage entries become NA", {
  zero_mat <- Matrix::sparseMatrix(
    i = integer(0), j = integer(0),
    dims = c(3L, 6L),
    dimnames = list(paste0("Cell", 1:3), paste0("jxn", 1:6))
  )
  result <- CalculatePSI(zero_mat, make_event_data(),
                          min_coverage = 1L, verbose = FALSE)
  expect_true(all(is.na(result)))
})

test_that("CalculatePSI (matrix): PSI = 1 when only inclusion reads present", {
  mat <- Matrix::sparseMatrix(
    i = c(1L, 1L), j = c(1L, 2L), x = c(5, 5),
    dims = c(1L, 6L),
    dimnames = list("Cell1", paste0("jxn", 1:6))
  )
  result <- CalculatePSI(mat, make_event_data()[1, ],
                          min_coverage = 1L, verbose = FALSE)
  expect_equal(result[1, 1], 1.0)
})

test_that("CalculatePSI (matrix): PSI = 0 when only exclusion reads present", {
  mat <- Matrix::sparseMatrix(
    i = 1L, j = 5L, x = 10,
    dims = c(1L, 6L),
    dimnames = list("Cell1", paste0("jxn", 1:6))
  )
  result <- CalculatePSI(mat, make_event_data()[1, ],
                          min_coverage = 1L, verbose = FALSE)
  expect_equal(result[1, 1], 0.0)
})

test_that("CalculatePSI (matrix): PSI = 0.5 with equal inclusion/exclusion", {
  mat <- Matrix::sparseMatrix(
    i = c(1L, 1L), j = c(1L, 5L), x = c(5, 5),
    dims = c(1L, 6L),
    dimnames = list("Cell1", paste0("jxn", 1:6))
  )
  result <- CalculatePSI(mat, make_event_data()[1, ],
                          min_coverage = 1L, verbose = FALSE)
  expect_equal(result[1, 1], 0.5)
})

test_that("CalculatePSI (matrix): min_coverage threshold is respected", {
  mat <- Matrix::sparseMatrix(
    i = c(1L, 1L, 2L, 2L), j = c(1L, 5L, 1L, 5L), x = c(2, 1, 5, 5),
    dims = c(2L, 6L),
    dimnames = list(c("Cell1", "Cell2"), paste0("jxn", 1:6))
  )
  result <- CalculatePSI(mat, make_event_data()[1, ],
                          min_coverage = 5L, verbose = FALSE)
  expect_true(is.na(result["Cell1", "SE_gene1_e2"]))
  expect_false(is.na(result["Cell2", "SE_gene1_e2"]))
})

test_that("CalculatePSI (matrix): handles missing junctions gracefully", {
  partial_mat <- Matrix::sparseMatrix(
    i = c(1L, 1L), j = c(1L, 2L), x = c(5, 5),
    dims = c(1L, 4L),
    dimnames = list("Cell1", paste0("jxn", 1:4))
  )
  result <- CalculatePSI(partial_mat, make_event_data()[1, ],
                          min_coverage = 1L, verbose = FALSE)
  expect_equal(result[1, 1], 1.0)
})

test_that("CalculatePSI (matrix): errors when events arg is missing", {
  jxn_mat <- make_junction_counts()
  expect_error(CalculatePSI(jxn_mat, verbose = FALSE))
})

# ---- PSI on a MatisseObject ------------------------------------------------

test_that("CalculatePSI (MatisseObject): populates 'psi' ChromatinAssay in Seurat", {
  skip_if_not_installed("Signac")
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  psi_assay <- GetSeurat(obj)[["psi"]]
  expect_false(is.null(psi_assay))
  expect_true(inherits(psi_assay, "ChromatinAssay"))
})

test_that("CalculatePSI (MatisseObject): GetPSI returns cells x events with correct dims", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  psi <- GetPSI(obj)
  expect_equal(nrow(psi), 10L)
  expect_equal(ncol(psi), 2L)
})

test_that("CalculatePSI (MatisseObject): inclusion and exclusion counts accessible", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  inc <- GetInclusionCounts(obj)
  exc <- GetExclusionCounts(obj)
  expect_false(is.null(inc))
  expect_false(is.null(exc))
  expect_equal(dim(inc), c(10L, 2L))
  expect_equal(dim(exc), c(10L, 2L))
})

test_that("CalculatePSI (MatisseObject): PSI values in [0,1] for covered entries", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, min_coverage = 1L, verbose = FALSE)
  psi <- as.matrix(GetPSI(obj))
  vals <- psi[!is.na(psi)]
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("CalculatePSI (MatisseObject): errors without junction_counts", {
  skip_if_not_installed("Seurat")
  seu <- make_seurat()
  obj <- CreateMatisseObject(seu, verbose = FALSE)
  expect_error(CalculatePSI(obj, verbose = FALSE),
               regexp = "junction_counts slot is NULL")
})

test_that("CalculatePSI (MatisseObject): errors without event_data", {
  skip_if_not_installed("Seurat")
  seu <- make_seurat()
  jxn <- make_junction_counts()
  obj <- CreateMatisseObject(seu, junction_counts = jxn, verbose = FALSE)
  expect_error(CalculatePSI(obj, verbose = FALSE),
               regexp = "No splice events defined")
})

# ---- SummarizePSI ----------------------------------------------------------

test_that("SummarizePSI: returns one row per event", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  s   <- SummarizePSI(obj)
  expect_equal(nrow(s), 2L)
  expect_true(all(c("event_id", "mean_psi", "median_psi", "sd_psi",
                    "n_cells_covered") %in% colnames(s)))
})

test_that("SummarizePSI: errors if 'psi' assay is absent", {
  obj <- make_matisse_object()
  expect_error(SummarizePSI(obj), regexp = "PSI matrix is empty")
})

# ---- CalculatePSI on transcript-based object (already has PSI) -------------

test_that("CalculatePSI: warns and returns unchanged when PSI already computed from transcripts", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Signac")
  obj <- make_matisse_from_transcripts()
  psi_before <- GetPSI(obj)
  expect_warning(
    obj2 <- CalculatePSI(obj, verbose = FALSE),
    regexp = "already has PSI computed"
  )
  # PSI should be unchanged
  expect_equal(as.matrix(GetPSI(obj2)), as.matrix(psi_before))
})

# ---- Internal PSI helpers --------------------------------------------------

test_that(".n_covered_per_cell: counts non-NA stored entries per row", {
  psi_sp <- Matrix::sparseMatrix(
    i = c(1, 1, 2), j = c(1, 3, 2), x = c(0.8, 0.2, 0.5),
    dims = c(3, 4), repr = "C")
  result <- Matisse:::.n_covered_per_cell(psi_sp)
  expect_equal(result, c(2L, 1L, 0L))
})

test_that(".psi_rowmeans_sparse: ignores NAs and absent entries", {
  psi_sp <- Matrix::sparseMatrix(
    i = c(1, 1, 2), j = c(1, 3, 2), x = c(0.8, 0.2, 0.5),
    dims = c(3, 4), repr = "C")
  result <- Matisse:::.psi_rowmeans_sparse(psi_sp)
  expect_equal(result[1], 0.5, tolerance = 1e-6)
  expect_equal(result[2], 0.5, tolerance = 1e-6)
  expect_true(is.na(result[3]))
})
