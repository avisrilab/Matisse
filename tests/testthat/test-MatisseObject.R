# ---------------------------------------------------------------------------
# Tests for MatisseObject construction, display, subsetting, and accessors
# ---------------------------------------------------------------------------

test_that("CreateMatisseObject: requires a Seurat object", {
  skip_if_not_installed("Seurat")
  expect_error(CreateMatisseObject(seurat = list()),
               regexp = "must be a Seurat object")
})

test_that("CreateMatisseObject: succeeds with minimal input", {
  obj <- make_matisse_object()
  expect_s4_class(obj, "MatisseObject")
})

test_that("CreateMatisseObject: mode is 'junction' when junction_counts supplied", {
  obj <- make_matisse_object()
  expect_equal(obj@mode, "junction")
})

test_that("CreateMatisseObject: mode is 'event' when ioe_files supplied", {
  obj <- make_matisse_from_transcripts()
  expect_equal(obj@mode, "event")
})

test_that("CreateMatisseObject: cell count matches Seurat object", {
  obj <- make_matisse_object()
  expect_equal(.n_cells(obj), 10L)
})

test_that("CreateMatisseObject: event_data is stored correctly", {
  obj <- make_matisse_object()
  ed  <- GetEventData(obj)
  expect_equal(nrow(ed), 2L)
  expect_true("event_id" %in% colnames(ed))
})

test_that("CreateMatisseObject: junction_data is stored correctly", {
  obj <- make_matisse_object()
  jd  <- GetJunctionData(obj)
  expect_equal(nrow(jd), 6L)
  expect_true(all(c("junction_id", "chr", "start", "end") %in% colnames(jd)))
})

test_that("CreateMatisseObject: junction counts stored as Assay5('junction')", {
  obj <- make_matisse_object()
  expect_false(is.null(GetSeurat(obj)[["junction"]]))
  expect_true(inherits(GetSeurat(obj)[["junction"]], "Assay5"))
})

test_that("CreateMatisseObject: GetJunctionCounts row names match cells", {
  obj <- make_matisse_object()
  jc  <- GetJunctionCounts(obj)
  expect_equal(rownames(jc), colnames(GetSeurat(obj)))
})

test_that("CreateMatisseObject: GetJunctionCounts returns cells x junctions", {
  obj <- make_matisse_object()
  jc  <- GetJunctionCounts(obj)
  expect_equal(nrow(jc), 10L)
  expect_equal(ncol(jc), 6L)
})

test_that("CreateMatisseObject: .n_junctions returns junction count in junction mode", {
  obj <- make_matisse_object()
  expect_equal(.n_junctions(obj), 6L)
})

test_that("CreateMatisseObject: .n_junctions returns 0 in event mode", {
  obj <- make_matisse_from_transcripts()
  expect_equal(.n_junctions(obj), 0L)
})

test_that("CreateMatisseObject: rejects event_data missing required columns", {
  skip_if_not_installed("Seurat")
  seu <- make_seurat()
  bad_events <- data.frame(event_id = "e1", gene_id = "g1")
  expect_error(
    CreateMatisseObject(seu, event_data = bad_events),
    regexp = "missing required columns"
  )
})

test_that("CreateMatisseObject: accepts transcript_counts and creates 'transcript' assay", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  obj    <- CreateMatisseObject(seu, transcript_counts = tx_mat, verbose = FALSE)
  expect_false(is.null(GetSeurat(obj)[["transcript"]]))
  tx <- GetTranscriptCounts(obj)
  expect_equal(ncol(tx), 10L)
})

test_that("show method: produces output without error", {
  obj <- make_matisse_object()
  expect_output(show(obj), regexp = "MatisseObject")
})

test_that("show method: junction-based mode label in show output", {
  obj <- make_matisse_object()
  expect_output(show(obj), regexp = "junction-based")
})

test_that("show method: event-based mode label in show output", {
  obj <- make_matisse_from_transcripts()
  expect_output(show(obj), regexp = "event-based")
})

test_that("show method: PSI coverage line present and not NA after CalculatePSI", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  out      <- capture.output(show(obj))
  cov_line <- grep("PSI coverage", out, value = TRUE)
  expect_length(cov_line, 1L)
  expect_false(grepl(":\\s+NA\\s+%", cov_line))
  expect_true(grepl("[0-9]", cov_line))
})

test_that("show method: lists 'psi' in Assays after CalculatePSI", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  out <- capture.output(show(obj))
  expect_true(any(grepl("psi", out)))
})

test_that("dim: returns c(n_cells, n_events) after PSI calculation", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  d   <- dim(obj)
  expect_equal(d[1], 10L)
  expect_equal(d[2], 2L)
})

test_that("dim: returns c(n_cells, 0) before PSI calculation", {
  obj <- make_matisse_object()
  expect_equal(dim(obj), c(10L, 0L))
})

test_that("subsetting [: reduces cell count", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  sub <- obj[paste0("Cell", 1:5), ]
  expect_equal(.n_cells(sub), 5L)
})

test_that("subsetting [: reduces event count in 'psi' assay", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  sub <- obj[, "SE-gene1-e2"]
  expect_equal(.n_events(sub), 1L)
})

test_that("subsetting [: PSI still accessible after event subset", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  sub <- obj[, "SE-gene1-e2"]
  psi <- GetPSI(sub)
  expect_equal(ncol(psi), 1L)
  expect_equal(colnames(psi), "SE-gene1-e2")
})

test_that("subsetting [: junction counts also subsetted", {
  obj <- make_matisse_object()
  sub <- obj[paste0("Cell", 1:5), ]
  jc  <- GetJunctionCounts(sub)
  expect_equal(nrow(jc), 5L)
  expect_equal(ncol(jc), 6L)
})

test_that("subsetting [: errors on unknown cell barcode", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  expect_error(obj["FAKE_CELL", ], regexp = "not found")
})

test_that("GetSeurat: returns an embedded Seurat object", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  expect_s4_class(GetSeurat(obj), "Seurat")
})

test_that("GetPSI: returns NULL before CalculatePSI", {
  obj <- make_matisse_object()
  expect_null(GetPSI(obj))
})

test_that("GetPSI: returns cells x events sparse matrix after CalculatePSI", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  psi <- GetPSI(obj)
  expect_true(inherits(psi, "Matrix"))
  expect_equal(nrow(psi), 10L)
  expect_equal(ncol(psi), 2L)
})

test_that("GetPSI: row names are cell barcodes, column names are event IDs", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  psi <- GetPSI(obj)
  expect_equal(rownames(psi), colnames(GetSeurat(obj)))
  expect_equal(colnames(psi), GetEventData(obj)$event_id)
})

test_that("GetInclusionCounts: returns cells x events matrix after CalculatePSI", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  inc <- GetInclusionCounts(obj)
  expect_false(is.null(inc))
  expect_equal(dim(inc), dim(GetPSI(obj)))
})

test_that("GetExclusionCounts: returns cells x events matrix after CalculatePSI", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  exc <- GetExclusionCounts(obj)
  expect_false(is.null(exc))
  expect_equal(dim(exc), dim(GetPSI(obj)))
})

test_that("PSI stored as 'psi' Assay5 in Seurat object", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  seu <- GetSeurat(obj)
  psi_assay <- seu[["psi"]]
  expect_false(is.null(psi_assay))
  expect_true(inherits(psi_assay, "Assay5"))
})

test_that("MatisseMeta: returns seurat@meta.data", {
  obj  <- make_matisse_object()
  meta <- MatisseMeta(obj)
  expect_true(is.data.frame(meta))
  expect_equal(nrow(meta), 10L)
  expect_true("orig.ident" %in% colnames(meta))
})

test_that("MatisseMeta<-: adds new column to seurat meta.data", {
  obj  <- make_matisse_object()
  cells <- colnames(GetSeurat(obj))
  new_df <- data.frame(my_col = 1:10, row.names = cells)
  MatisseMeta(obj) <- new_df
  expect_true("my_col" %in% colnames(MatisseMeta(obj)))
})

test_that("AddIsoformMetadata: adds new columns to seurat meta.data", {
  obj    <- make_matisse_object()
  cells  <- colnames(GetSeurat(obj))
  new_df <- data.frame(my_new_col = as.numeric(1:10), row.names = cells)
  obj    <- AddIsoformMetadata(obj, new_df)
  expect_equal(nrow(MatisseMeta(obj)), 10L)
  expect_true("my_new_col" %in% colnames(MatisseMeta(obj)))
})

test_that("$ operator: returns Seurat meta.data column", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  expect_no_error(obj$orig.ident)
})

test_that(".DollarNames: lists Seurat metadata columns for tab-completion", {
  obj <- make_matisse_with_umap()   # has cell_type and seurat_clusters
  nms <- .DollarNames(obj, pattern = "")
  expect_true("cell_type" %in% nms)
  expect_true("seurat_clusters" %in% nms)
})

test_that(".DollarNames: pattern filters returned names", {
  obj <- make_matisse_with_umap()
  nms <- .DollarNames(obj, pattern = "cell")
  expect_true(all(grepl("cell", nms)))
})

test_that("$ operator: returns a Seurat function as a forwarding closure", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  fn  <- obj$NormalizeData
  expect_true(is.function(fn))
})

test_that("$ operator: calling returned closure runs on embedded Seurat and returns MatisseObject", {
  skip_if_not_installed("Seurat")
  obj    <- make_matisse_object()
  result <- obj$NormalizeData()
  expect_s4_class(result, "MatisseObject")
})

test_that("SetPSI: updates the data layer in the 'psi' assay", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  old_psi <- GetPSI(obj)
  new_psi <- old_psi * 0.5
  obj     <- SetPSI(obj, new_psi)
  updated <- GetPSI(obj)
  expect_equal(as.matrix(updated), as.matrix(new_psi), tolerance = 1e-6)
})

test_that("GetJunctionCounts: returns NULL for event-mode objects", {
  obj <- make_matisse_from_transcripts()
  expect_null(GetJunctionCounts(obj))
})

test_that("validity: passes for a valid object", {
  skip_if_not_installed("Seurat")
  seu2 <- make_seurat(n_cells = 5L)
  obj  <- methods::new("MatisseObject",
                        seurat     = seu2,
                        event_data = make_event_data())
  expect_no_error(methods::validObject(obj))
})
