# ---------------------------------------------------------------------------
# Tests for PSI calculation (CalculatePSI, SummarizePSI, internal helpers)
# ---------------------------------------------------------------------------

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
    dimnames = list("Cell1", .junction_ids_fixture)
  )
  result <- CalculatePSI(mat, make_event_data()[1, ],
                          min_coverage = 1L, verbose = FALSE)
  expect_equal(result[1, 1], 1.0)
})

test_that("CalculatePSI (matrix): PSI = 0 when only exclusion reads present", {
  mat <- Matrix::sparseMatrix(
    i = 1L, j = 5L, x = 10,
    dims = c(1L, 6L),
    dimnames = list("Cell1", .junction_ids_fixture)
  )
  result <- CalculatePSI(mat, make_event_data()[1, ],
                          min_coverage = 1L, verbose = FALSE)
  expect_equal(result[1, 1], 0.0)
})

test_that("CalculatePSI (matrix): PSI = 0.5 with symmetric SE reads (per-junction normalization)", {
  # PSI is normalized per junction (SUPPA2 / rMATS convention). For an SE
  # event the fixture has 2 inclusion junctions and 1 exclusion junction.
  # "Equal inclusion/exclusion" therefore means equal *per-junction* read
  # counts: inc_up = inc_dn = exc. With 5 reads at each, the molecular
  # fraction is 0.5 (each inclusion molecule produces reads at both
  # inclusion junctions, each exclusion molecule produces reads at one).
  mat <- Matrix::sparseMatrix(
    i = c(1L, 1L, 1L), j = c(1L, 2L, 5L), x = c(5, 5, 5),
    dims = c(1L, 6L),
    dimnames = list("Cell1", .junction_ids_fixture)
  )
  result <- CalculatePSI(mat, make_event_data()[1, ],
                          min_coverage = 1L, verbose = FALSE)
  expect_equal(result[1, 1], 0.5)
})

test_that("CalculatePSI (matrix): SE event with 2 inc + 1 exc junctions normalizes per junction", {
  # Regression for the bias the old raw-sum formula introduced.
  #
  # Biological setup: an SE (skipped exon) event has 2 inclusion junctions
  # (J_up: E1->cassette, J_dn: cassette->E3) and 1 exclusion junction
  # (J_skip: E1->E3). One inclusion-isoform molecule produces reads at
  # BOTH inclusion junctions over time; one exclusion-isoform molecule
  # produces reads at ONE exclusion junction.
  #
  # If true PSI = 0.5 (equal numbers of inclusion and exclusion molecules)
  # and coverage is uniform, you observe roughly equal reads at every
  # junction: J_up = J_dn = J_skip = 50.
  #
  # The naive raw-sum formula (inc_up+inc_dn) / (inc_up+inc_dn+exc) gives
  # 100/150 = 0.667 — biased toward inclusion because the same molecule
  # is counted twice (once at each inclusion junction). Per-junction
  # normalization (mean-of-inc / (mean-of-inc + mean-of-exc)) gives 0.5,
  # which is the molecular fraction users expect.
  events <- data.frame(
    event_id           = "SE_test",
    gene_id            = "G1",
    chr                = "chr1",
    strand             = "+",
    event_type         = "SE",
    inclusion_features = "chr1-100-200-+;chr1-300-400-+",
    exclusion_features = "chr1-100-400-+",
    stringsAsFactors   = FALSE
  )
  mat <- Matrix::Matrix(
    matrix(c(50, 50, 50), nrow = 1, ncol = 3,
           dimnames = list("Cell1",
                           c("chr1-100-200-+", "chr1-300-400-+",
                             "chr1-100-400-+"))),
    sparse = TRUE
  )
  result <- CalculatePSI(mat, events, min_coverage = 1L, verbose = FALSE)
  expect_equal(result["Cell1", "SE_test"], 0.5, tolerance = 1e-9)
})

test_that("CalculatePSI (matrix): min_coverage threshold is respected", {
  mat <- Matrix::sparseMatrix(
    i = c(1L, 1L, 2L, 2L), j = c(1L, 5L, 1L, 5L), x = c(2, 1, 5, 5),
    dims = c(2L, 6L),
    dimnames = list(c("Cell1", "Cell2"), .junction_ids_fixture)
  )
  result <- CalculatePSI(mat, make_event_data()[1, ],
                          min_coverage = 5L, verbose = FALSE)
  expect_true(is.na(result["Cell1", "SE-gene1-e2"]))
  expect_false(is.na(result["Cell2", "SE-gene1-e2"]))
})

test_that("CalculatePSI (matrix): handles missing junctions gracefully", {
  partial_mat <- Matrix::sparseMatrix(
    i = c(1L, 1L), j = c(1L, 2L), x = c(5, 5),
    dims = c(1L, 4L),
    dimnames = list("Cell1", .junction_ids_fixture[1:4])
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

test_that("CalculatePSI (MatisseObject): populates 'psi' Assay5 in Seurat", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  psi_assay <- GetSeurat(obj)[["psi"]]
  expect_false(is.null(psi_assay))
  expect_true(inherits(psi_assay, "Assay5"))
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
  inc <- Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "counts"))
  exc <- Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "exclusion"))
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

test_that("CalculatePSI (MatisseObject): errors without junction assay", {
  skip_if_not_installed("Seurat")
  seu <- make_seurat()
  obj <- CreateMatisseObject(seu, verbose = FALSE)
  expect_error(CalculatePSI(obj, verbose = FALSE),
               regexp = "No junction assay found")
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

test_that("SummarizePSI: no NaN in output (uncovered events yield NA)", {
  # Build an object where min_coverage is very high so no event is covered
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, min_coverage = 9999L, verbose = FALSE)
  s   <- SummarizePSI(obj)
  expect_false(any(is.nan(s$mean_psi)))
  expect_false(any(is.nan(s$median_psi)))
  expect_false(any(is.nan(s$sd_psi)))
  expect_true(all(is.na(s$mean_psi)))
})

# ---- CalculatePSI can be called repeatedly (idempotent) --------------------

test_that("CalculatePSI: can be called again on an object that already has PSI", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_from_transcripts()
  expect_no_error(CalculatePSI(obj, verbose = FALSE))
})

# ---- Internal PSI helpers --------------------------------------------------

test_that(".n_covered_per_cell: counts non-NA stored entries per row", {
  psi_sp <- Matrix::sparseMatrix(
    i = c(1, 1, 2), j = c(1, 3, 2), x = c(0.8, 0.2, 0.5),
    dims = c(3, 4), repr = "C")
  result <- Matisse:::.n_covered_per_cell(psi_sp)
  expect_equal(result, c(2L, 1L, 0L))
})

test_that(".build_indicator_matrix: matches across SeuratObject's `_` -> `-` sanitisation", {
  # Regression: SeuratObject sanitises feature names by replacing
  # underscores with dashes when building an Assay5. IOE inclusion /
  # exclusion feature columns keep their underscores, so a raw
  # match() against the post-sanitisation rownames missed every entry
  # and silently zeroed all per-event read counts in transcript mode
  # for any IOE input with underscored transcript IDs (RefSeq style).
  # .build_indicator_matrix now sanitises both sides before matching.
  id_lists    <- list(c("tx_inc_a", "tx_inc_b"), c("tx_inc_c"))
  id_universe <- c("tx-inc-a", "tx-inc-b", "tx-inc-c", "tx-exc-a")
  A <- Matisse:::.build_indicator_matrix(id_lists, id_universe)
  expect_equal(as.numeric(A[, 1]), c(1, 1, 0, 0))
  expect_equal(as.numeric(A[, 2]), c(0, 0, 1, 0))
})

test_that("CalculatePSI (transcript mode): underscored transcript IDs aggregate correctly (regression)", {
  # End-to-end check that the underscore -> dash sanitisation fix lets a
  # full IOE input with underscored transcript IDs flow through
  # CreateMatisseObject + CalculatePSI without zeroing the inc/exc
  # layers.
  obj <- make_matisse_long_read()  # SE event with tx_inc_a etc.
  psi_assay <- GetSeurat(obj)[["psi"]]
  inc <- as.numeric(SeuratObject::GetAssayData(
    psi_assay, layer = "counts")["SE:chr1:1201-2999:3201-4999:+", ])
  exc <- as.numeric(SeuratObject::GetAssayData(
    psi_assay, layer = "exclusion")["SE:chr1:1201-2999:3201-4999:+", ])
  expect_gt(sum(inc), 0)
  expect_gt(sum(exc), 0)
})

