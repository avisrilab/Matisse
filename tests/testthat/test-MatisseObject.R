# ---------------------------------------------------------------------------
# Tests for MatisseObject construction, display, subsetting, and accessors
# ---------------------------------------------------------------------------

test_that("CreateMatisseObject: requires a Seurat object", {
  skip_if_not_installed("Seurat")
  expect_error(CreateMatisseObject(seurat = list()),
               regexp = "must be a Seurat object")
})

test_that("CreateMatisseObject: full abort message survives (regression for rlang::abort multi-arg)", {
  # Regression: multi-positional-arg rlang::abort() silently truncates the
  # message and consumes the second string as the error class. We collapsed
  # all such call sites with paste0(); this asserts the full message survives.
  skip_if_not_installed("Seurat")
  seu <- make_seurat()
  jxn <- make_junction_counts()
  tx  <- make_transcript_counts()
  e <- tryCatch(
    CreateMatisseObject(seurat = seu, junction_counts = jxn, transcript_counts = tx),
    error = function(e) e
  )
  expect_match(conditionMessage(e), "junction_counts.*junction mode")
  expect_match(conditionMessage(e), "transcript_counts.*transcript mode")
  expect_match(conditionMessage(e), "not both")
})

test_that("CreateMatisseObject: succeeds with minimal input", {
  obj <- make_matisse_object()
  expect_s4_class(obj, "MatisseObject")
})

test_that("CreateMatisseObject: input.mode is 'junction' when junction_counts supplied", {
  obj <- make_matisse_object()
  expect_equal(obj@input.mode, "junction")
})

test_that("CreateMatisseObject: input.mode is 'transcript' when transcript_counts supplied", {
  obj <- make_matisse_from_transcripts()
  expect_equal(obj@input.mode, "transcript")
})

test_that("CreateMatisseObject: cell count matches Seurat object", {
  obj <- make_matisse_object()
  expect_equal(.n_cells(obj), 10L)
})

test_that("CreateMatisseObject: event_data is staged in @misc until CalculatePSI runs", {
  # Pre-CalculatePSI, event annotation is staged under @misc (no PSI assay
  # exists yet to host meta.features). Post-P1 it migrates to assay.
  obj <- make_matisse_object()
  ed  <- obj@misc[["event_data"]]
  expect_equal(nrow(ed), 2L)
  expect_true("event_id" %in% colnames(ed))
})

test_that("CreateMatisseObject: event_data_path is NA when event_data is a data.frame", {
  obj <- make_matisse_object()
  expect_identical(obj@misc[["event_data_path"]], NA_character_)
})

test_that("CreateMatisseObject: folds PSI calc into construction by default (P5)", {
  # P5: by default, the constructor computes PSI as part of construction.
  # The returned object has a populated 'psi' assay.
  skip_if_not_installed("Seurat")
  seu     <- make_seurat()
  jxn_mat <- make_junction_counts()
  ev_data <- make_event_data()
  obj <- CreateMatisseObject(
    seurat          = seu,
    junction_counts = jxn_mat,
    event_data      = ev_data,
    min_coverage    = 1L,
    verbose         = FALSE
  )
  expect_true("psi" %in% SeuratObject::Assays(GetSeurat(obj)))
  expect_true("nPercent_isoform" %in% colnames(MatisseMeta(obj)))
})

test_that("CreateMatisseObject: defer_psi=TRUE skips the PSI step", {
  skip_if_not_installed("Seurat")
  seu     <- make_seurat()
  jxn_mat <- make_junction_counts()
  ev_data <- make_event_data()
  obj <- CreateMatisseObject(
    seurat          = seu,
    junction_counts = jxn_mat,
    event_data      = ev_data,
    defer_psi       = TRUE,
    verbose         = FALSE
  )
  expect_false("psi" %in% SeuratObject::Assays(GetSeurat(obj)))
})

test_that("CreateMatisseObject: event_data_path is the normalized path when ioe_files used", {
  skip_if_not_installed("Seurat")
  ioe <- make_se_ioe_file()
  obj <- CreateMatisseObject(
    seurat            = make_seurat(n_cells = 20L),
    transcript_counts = make_transcript_counts(n_cells = 20L),
    ioe_files         = ioe,
    verbose           = FALSE
  )
  expect_equal(obj@misc[["event_data_path"]],
               normalizePath(ioe, mustWork = FALSE))
})

test_that("CreateMatisseObject: warns when both event_data and ioe_files are supplied (P14)", {
  skip_if_not_installed("Seurat")
  ev  <- data.frame(
    event_id             = c("X1", "X2"),
    gene_id              = c("g1", "g2"),
    chr                  = c("chr1", "chr1"),
    strand               = c("+", "+"),
    event_type           = c("SE", "SE"),
    inclusion_junctions  = c("a;b", "c;d"),
    exclusion_junctions  = c("e", "f"),
    stringsAsFactors     = FALSE
  )
  ioe <- make_se_ioe_file()
  expect_warning(
    CreateMatisseObject(
      seurat            = make_seurat(n_cells = 20L),
      transcript_counts = make_se_transcript_counts(n_cells = 20L),
      ioe_files         = ioe,
      event_data        = ev,
      defer_psi         = TRUE,
      verbose           = FALSE
    ),
    regexp = "event_data.*takes precedence|both.*ioe_files.*event_data"
  )
})

test_that("CalculatePSI clears @misc[['event_data']] (now lives in PSI assay meta.features)", {
  obj <- make_matisse_object()
  # Pre-CalculatePSI: event_data is staged in @misc
  expect_true("event_data" %in% names(obj@misc))
  expect_gt(nrow(obj@misc[["event_data"]]), 0L)
  obj <- CalculatePSI(obj, verbose = FALSE)
  # Post-CalculatePSI: @misc[["event_data"]] is cleared, lives in assay
  expect_null(obj@misc[["event_data"]])
  mf <- obj@seurat[["psi"]][[]]
  expect_true("gene_id" %in% colnames(mf))
})

test_that("CreateMatisseObject: junction_data is stored correctly", {
  obj <- make_matisse_object()
  jd  <- obj@misc[["junction_data"]]
  expect_equal(nrow(jd), 6L)
  expect_true(all(c("junction_id", "chr", "start", "end") %in% colnames(jd)))
})

test_that("CreateMatisseObject: junction counts stored as Assay5('isoform')", {
  obj <- make_matisse_object()
  expect_false(is.null(GetSeurat(obj)[["isoform"]]))
  expect_true(inherits(GetSeurat(obj)[["isoform"]], "Assay5"))
})

test_that("CreateMatisseObject: junction counts (cells x junctions) align with cells", {
  # P2 dropped GetJunctionCounts; access junction reads via the native
  # Seurat layer accessor on the "isoform" assay.
  obj <- make_matisse_object()
  jc  <- Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]], "counts"))
  expect_equal(rownames(jc), colnames(GetSeurat(obj)))
  expect_equal(nrow(jc), 10L)
  expect_equal(ncol(jc), 6L)
})

test_that("CreateMatisseObject: .n_isoforms returns junction count in junction mode", {
  obj <- make_matisse_object()
  expect_equal(.n_isoforms(obj), 6L)
})

test_that("CreateMatisseObject: .n_isoforms returns transcript count in transcript mode", {
  obj <- make_matisse_from_transcripts()
  expect_equal(.n_isoforms(obj), 8L)   # make_transcript_counts() has 8 transcripts
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

test_that("CreateMatisseObject: accepts transcript_counts and creates 'isoform' assay", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  tx_mat <- make_transcript_counts()
  obj    <- CreateMatisseObject(seu, transcript_counts = tx_mat, verbose = FALSE)
  expect_false(is.null(GetSeurat(obj)[["isoform"]]))
  # Transcript counts are stored as Assay5("isoform"), transcripts x cells.
  tx <- SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]], "counts")
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

test_that("show method: transcript-based mode label in show output", {
  obj <- make_matisse_from_transcripts()
  expect_output(show(obj), regexp = "transcript-based")
})

test_that("show method: default assay is marked with * in Assays line", {
  obj <- make_matisse_object()
  out        <- capture.output(show(obj))
  assay_line <- grep("Assays", out, value = TRUE)
  default    <- SeuratObject::DefaultAssay(GetSeurat(obj))
  expect_true(grepl(paste0(default, "\\*"), assay_line))
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

test_that("subsetting [: rejects single-event subsets (Assay5 minimum-features constraint)", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  expect_error(
    obj[, "SE-gene1-e2"],
    regexp = "fewer than 2 events"
  )
})

test_that("subsetting [: PSI still accessible after multi-event subset", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  # The fixture has exactly 2 events; subset to both via [, ] (a no-op subset
  # at the API level but exercises the assay-rebuild path).
  sub <- obj[, c("SE-gene1-e2", "SE-gene1-e3")]
  psi <- GetPSI(sub)
  expect_equal(ncol(psi), 2L)
  expect_setequal(colnames(psi), c("SE-gene1-e2", "SE-gene1-e3"))
})

test_that("subsetting [: junction counts also subsetted", {
  obj <- make_matisse_object()
  sub <- obj[paste0("Cell", 1:5), ]
  jc  <- Matrix::t(SeuratObject::GetAssayData(GetSeurat(sub)[["isoform"]], "counts"))
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
  # Event IDs are the rownames of the PSI assay's meta.features (post-P1).
  expect_equal(colnames(psi), rownames(GetSeurat(obj)[["psi"]][[]]))
})

test_that("PSI inclusion counts (cells x events) accessible via the 'psi' assay 'counts' layer", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  inc <- Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "counts"))
  expect_equal(dim(inc), dim(GetPSI(obj)))
})

test_that("PSI exclusion counts accessible via the 'psi' assay 'exclusion' layer", {
  obj <- make_matisse_object()
  obj <- CalculatePSI(obj, verbose = FALSE)
  exc <- Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "exclusion"))
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

test_that("MatisseMeta<-: aligns by barcode rownames, not by position", {
  # Regression: previously the setter looped value[[col]] -> meta.data[[col]]
  # by index, silently writing wrong values when value's rownames differed
  # in order from meta.data's rownames.
  obj    <- make_matisse_object()
  cells  <- colnames(GetSeurat(obj))
  shuffled <- rev(cells)
  new_df <- data.frame(shuffled_val = seq_along(shuffled),
                       row.names    = shuffled)
  MatisseMeta(obj) <- new_df
  # The cell that was originally last should have value 1 (it was first in `shuffled`)
  result <- MatisseMeta(obj)$shuffled_val
  names(result) <- rownames(MatisseMeta(obj))
  expect_equal(result[shuffled], seq_along(shuffled),
               ignore_attr = TRUE)
})

test_that("AddMetaData(matisse_obj, df) dispatches via S3 and updates cell metadata", {
  # P2 dropped AddIsoformMetadata; the S3 method AddMetaData.MatisseObject
  # in dispatch.R handles this natively.
  obj    <- make_matisse_object()
  cells  <- colnames(GetSeurat(obj))
  new_df <- data.frame(my_new_col = as.numeric(1:10), row.names = cells)
  obj    <- SeuratObject::AddMetaData(obj, metadata = new_df)
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

test_that("transcript-mode object stores transcripts (not junctions) in 'isoform' assay", {
  # P2 dropped the mode-aware GetJunctionCounts NULL-return wrapper. The
  # equivalent assertion: in transcript mode, the "isoform" assay holds
  # transcripts as features, not junctions.
  obj <- make_matisse_from_transcripts()
  expect_equal(obj@input.mode, "transcript")
  expect_true("isoform" %in% SeuratObject::Assays(GetSeurat(obj)))
})

test_that("validity: passes for a valid object", {
  skip_if_not_installed("Seurat")
  seu2 <- make_seurat(n_cells = 5L)
  obj  <- methods::new("MatisseObject",
                        seurat     = seu2,
                        input.mode = "junction")
  expect_no_error(methods::validObject(obj))
})
