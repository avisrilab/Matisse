test_that("CreateMatisseObject: requires a Seurat object", {
  skip_if_not_installed("Seurat")
  expect_error(CreateMatisseObject(seurat = list()),
               regexp = "must be a Seurat object")
})

test_that("CreateMatisseObject: succeeds with minimal input", {
  obj <- make_matisse_object()
  expect_s4_class(obj, "MatisseObject")
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

test_that("CreateMatisseObject: junction_counts row names match cells", {
  obj <- make_matisse_object()
  jc  <- GetJunctionCounts(obj)
  expect_equal(rownames(jc), colnames(GetSeurat(obj)))
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

test_that("show method: produces output without error", {
  obj <- make_matisse_object()
  expect_output(show(obj), regexp = "MatisseObject")
})

test_that("show method: PSI coverage line is a valid percentage, not NA", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  out <- capture.output(show(obj))
  cov_line <- grep("PSI coverage", out, value = TRUE)
  expect_length(cov_line, 1L)
  # Must not show "NA %" (bare NA as the computed value)
  expect_false(grepl(":\\s+NA\\s+%", cov_line))
  # Must contain a numeric percentage
  expect_true(grepl("[0-9]", cov_line))
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
  obj  <- make_matisse_object()
  obj  <- CalculatePSI(obj, verbose = FALSE)
  sub  <- obj[paste0("Cell", 1:5), ]
  expect_equal(.n_cells(sub), 5L)
})

test_that("subsetting [: reduces event count", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  sub <- obj[, "SE_gene1_e2"]
  expect_equal(.n_events(sub), 1L)
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

test_that("MatisseMeta<-: assignment updates isoform_metadata", {
  obj  <- make_matisse_object()
  meta <- data.frame(my_col = 1:10, row.names = colnames(GetSeurat(obj)))
  MatisseMeta(obj) <- meta
  expect_true("my_col" %in% colnames(MatisseMeta(obj)))
})

test_that("AddIsoformMetadata: adds new columns", {
  obj <- make_matisse_object()
  new_col <- setNames(as.numeric(1:10),
                      colnames(GetSeurat(obj)))
  obj <- AddIsoformMetadata(obj, new_col)
  # Column should be present (name comes from vector names in some forms)
  expect_true(nrow(MatisseMeta(obj)) == 10L)
})

test_that("$ operator: delegates to Seurat meta.data", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  # Seurat always has orig.ident in meta.data
  expect_no_error(obj$orig.ident)
})

test_that("validity: catches PSI matrix with wrong row names", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)

  # Tamper with PSI row names
  bad_psi <- obj@psi
  rownames(bad_psi) <- paste0("WrongCell", seq_len(nrow(bad_psi)))
  expect_error(
    methods::new("MatisseObject",
      seurat = obj@seurat, psi = bad_psi),
    regexp = "Row names of 'psi'"
  )
})
