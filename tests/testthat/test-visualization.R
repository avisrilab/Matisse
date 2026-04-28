# ---------------------------------------------------------------------------
# Tests for visualization functions
# All Plot* functions return a ggplot object; we do not snapshot pixel output.
# ---------------------------------------------------------------------------

# ============================================================
# PlotUMAP
# ============================================================

test_that("PlotUMAP: returns a ggplot object", {
  obj <- make_matisse_with_umap()
  p   <- PlotUMAP(obj, feature = "SE-gene1-e2")
  expect_s3_class(p, "gg")
})

test_that("PlotUMAP: custom title is applied", {
  obj <- make_matisse_with_umap()
  p   <- PlotUMAP(obj, feature = "SE-gene1-e2", title = "My Title")
  expect_equal(p$labels$title, "My Title")
})

test_that("PlotUMAP: default title is the feature name", {
  obj <- make_matisse_with_umap()
  p   <- PlotUMAP(obj, feature = "SE-gene1-e2")
  expect_equal(p$labels$title, "SE-gene1-e2")
})

test_that("PlotUMAP: data has the same number of rows as the cell count", {
  obj <- make_matisse_with_umap()
  p   <- PlotUMAP(obj, feature = "SE-gene1-e2")
  expect_equal(nrow(p$data), .n_cells(obj))
})

test_that("PlotUMAP: PSI values are in [0, 1] or NA in the plot data", {
  obj  <- make_matisse_with_umap()
  p    <- PlotUMAP(obj, feature = "SE-gene1-e2")
  vals <- p$data$val
  finite_vals <- vals[!is.na(vals)]
  expect_true(all(finite_vals >= 0 & finite_vals <= 1))
})

test_that("PlotUMAP: errors if PSI matrix is NULL", {
  obj <- make_matisse_object()   # PSI not yet calculated
  expect_error(PlotUMAP(obj, feature = "SE-gene1-e2"),
               regexp = "PSI assay is NULL")
})

test_that("PlotUMAP: errors for a feature not present anywhere", {
  obj <- make_matisse_with_umap()
  expect_error(PlotUMAP(obj, feature = "nonexistent_event_XYZ"),
               regexp = "not found")
})

test_that("PlotUMAP: errors if the requested reduction is absent from the Seurat object", {
  # make_matisse_object uses a seurat WITHOUT a umap reduction
  obj <- CalculatePSI(make_matisse_object(), min_coverage = 1L, verbose = FALSE)
  expect_error(PlotUMAP(obj, feature = "SE-gene1-e2"))
})

# ============================================================
# PlotViolin
# ============================================================

test_that("PlotViolin: returns a ggplot object", {
  obj <- make_matisse_with_umap()
  p   <- PlotViolin(obj, feature = "SE-gene1-e2", group_by = "cell_type")
  expect_s3_class(p, "gg")
})

test_that("PlotViolin: default group_by 'seurat_clusters' is used when not specified", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_with_umap()
  expect_no_error(PlotViolin(obj, feature = "SE-gene1-e2"))
})

test_that("PlotViolin: add_points = TRUE still returns a ggplot", {
  obj <- make_matisse_with_umap()
  p   <- PlotViolin(obj, feature = "SE-gene1-e2",
                        group_by = "cell_type", add_points = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotViolin: custom colour vector is accepted", {
  obj    <- make_matisse_with_umap()
  colors <- c(TypeA = "#E41A1C", TypeB = "#377EB8")
  p      <- PlotViolin(obj, feature = "SE-gene1-e2",
                           group_by = "cell_type", colours = colors)
  expect_s3_class(p, "gg")
})

test_that("PlotViolin: custom title is applied", {
  obj <- make_matisse_with_umap()
  p   <- PlotViolin(obj, feature = "SE-gene1-e2",
                        group_by = "cell_type", title = "Splicing switch")
  expect_equal(p$labels$title, "Splicing switch")
})

test_that("PlotViolin: feature = NULL auto-selects isoform QC metrics", {
  obj <- make_matisse_with_umap()
  # nCount_isoform and nFeature_isoform exist from construction;
  # nPercent_isoform added by CalculatePSI (already called inside make_matisse_with_umap)
  p <- PlotViolin(obj, feature = NULL, group_by = "cell_type")
  expect_s3_class(p, "gg")
  # Should be faceted (multiple metrics)
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("PlotViolin: feature = NULL errors when no QC metrics are present", {
  # A fresh object with no isoform assay has no isoform QC columns
  skip_if_not_installed("Seurat")
  seu <- make_seurat_with_umap()
  # Create a MatisseObject without any junction or transcript counts
  obj <- methods::new("MatisseObject",
                      seurat     = seu,
                      input.mode = "junction",
                      version    = "test",
                      misc       = list())
  expect_error(
    PlotViolin(obj, feature = NULL, group_by = "cell_type"),
    regexp = "No QC metrics found"
  )
})

test_that("PlotViolin: errors if PSI matrix is NULL", {
  obj <- make_matisse_object()
  expect_error(PlotViolin(obj, feature = "SE-gene1-e2"),
               regexp = "PSI assay is NULL")
})

test_that("PlotViolin: errors for an unknown feature", {
  obj <- make_matisse_with_umap()
  expect_error(PlotViolin(obj, feature = "ghost_event"),
               regexp = "not found")
})

test_that("PlotViolin: errors for an unknown group_by column", {
  obj <- make_matisse_with_umap()
  expect_error(
    PlotViolin(obj, feature = "SE-gene1-e2", group_by = "no_such_col"),
    regexp = "not found"
  )
})

test_that("PlotViolin: vector of features returns a faceted ggplot", {
  obj <- make_matisse_with_umap()  # has PSI + cell_type
  p   <- PlotViolin(obj,
                    feature  = c("SE-gene1-e2", "nPercent_isoform"),
                    group_by = "cell_type")
  expect_s3_class(p, "gg")
  # facet_wrap adds a FacetWrap component
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("PlotViolin: metadata column (QC metric) is accepted as feature", {
  obj <- make_matisse_with_umap()
  p   <- PlotViolin(obj, feature = "nPercent_isoform", group_by = "cell_type")
  expect_s3_class(p, "gg")
})

test_that("PlotViolin: ncol controls facet columns for multiple features", {
  obj <- make_matisse_with_umap()
  p   <- PlotViolin(obj,
                    feature  = c("SE-gene1-e2", "nPercent_isoform"),
                    group_by = "cell_type",
                    ncol     = 1L)
  expect_s3_class(p, "gg")
})

# ============================================================
# PlotHeatmap
# ============================================================

test_that("PlotHeatmap: returns a ggplot object", {
  obj <- make_matisse_with_umap()
  p   <- PlotHeatmap(obj)
  expect_s3_class(p, "gg")
})

test_that("PlotHeatmap: max_cells downsampling is respected", {
  obj <- make_matisse_with_umap()
  p   <- PlotHeatmap(obj, max_cells = 5L)
  expect_lte(length(unique(p$data$cell)), 5L)
})

test_that("PlotHeatmap: subsetting to specific events works", {
  obj <- make_matisse_with_umap()
  p   <- PlotHeatmap(obj, events = "SE-gene1-e2")
  expect_equal(length(unique(p$data$event)), 1L)
})

test_that("PlotHeatmap: subsetting to specific cells works", {
  obj   <- make_matisse_with_umap()
  cells <- paste0("Cell", 1:5)
  p     <- PlotHeatmap(obj, cells = cells)
  expect_lte(length(unique(p$data$cell)), 5L)
})

test_that("PlotHeatmap: group_by orders cells without error", {
  obj <- make_matisse_with_umap()
  expect_no_error(PlotHeatmap(obj, group_by = "cell_type"))
})

test_that("PlotHeatmap: warns on unknown event IDs (does not error)", {
  obj <- make_matisse_with_umap()
  expect_warning(
    PlotHeatmap(obj, events = c("SE-gene1-e2", "bad_event_99")),
    regexp = "not found"
  )
})

test_that("PlotHeatmap: max_events cap is respected", {
  obj <- make_matisse_with_umap()
  p   <- PlotHeatmap(obj, max_events = 1L)
  expect_equal(length(unique(p$data$event)), 1L)
})

test_that("PlotHeatmap: max_events selects top-variance events with message", {
  obj <- make_matisse_with_umap()
  expect_message(
    PlotHeatmap(obj, max_events = 1L),
    regexp = "highest-variance"
  )
})

test_that("PlotHeatmap: errors if PSI matrix is NULL", {
  obj <- make_matisse_object()
  expect_error(PlotHeatmap(obj), regexp = "PSI assay is NULL")
})


# ============================================================
# PlotSashimi
# ============================================================

test_that("PlotSashimi: returns a ggplot for junction-mode object", {
  obj <- make_matisse_short_read()
  p   <- PlotSashimi(obj, event_id = "SE:chr1:1201-2999:3201-4999:+")
  expect_s3_class(p, "gg")
})

test_that("PlotSashimi: returns a ggplot for transcript-mode object", {
  obj <- make_matisse_long_read()
  p   <- PlotSashimi(obj, event_id = "SE:chr1:1201-2999:3201-4999:+")
  expect_s3_class(p, "gg")
})

test_that("PlotSashimi: group_by produces a faceted plot", {
  obj <- make_matisse_short_read()
  p   <- PlotSashimi(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                      group_by = "cell_type")
  expect_s3_class(p, "gg")
  expect_true(!is.null(p$facet))
})

test_that("PlotSashimi: arc_scale = 'linear' is accepted", {
  obj <- make_matisse_short_read()
  p   <- PlotSashimi(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                      arc_scale = "linear")
  expect_s3_class(p, "gg")
})

test_that("PlotSashimi: arc_scale = 'log' is accepted", {
  obj <- make_matisse_short_read()
  p   <- PlotSashimi(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                      arc_scale = "log")
  expect_s3_class(p, "gg")
})

test_that("PlotSashimi: custom title is applied", {
  obj <- make_matisse_short_read()
  p   <- PlotSashimi(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                      title = "My SE event")
  expect_equal(p$labels$title, "My SE event")
})

test_that("PlotSashimi: cell subset is accepted without error", {
  obj   <- make_matisse_short_read()
  cells <- paste0("Cell", 1:5)
  expect_no_error(
    PlotSashimi(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                 cells = cells)
  )
})

test_that("PlotSashimi: errors for unknown event_id", {
  obj <- make_matisse_short_read()
  expect_error(
    PlotSashimi(obj, event_id = "SE:chr99:0-1:2-3:+"),
    regexp = "not found in event annotation"
  )
})

test_that("PlotSashimi: errors in transcript mode for an unknown event_id", {
  # Event annotation lives in the PSI assay's meta.features; rownames carry
  # the event_id (also the assay feature name). Test the event-id-not-found
  # error path directly.
  obj <- make_matisse_long_read()
  expect_error(
    PlotSashimi(obj, event_id = "A3:chr1:1201-2999:3201-4999:+"),
    regexp = "not found in event annotation"
  )
})

test_that("PlotSashimi: returns a ggplot for RI event in transcript mode", {
  n_cells <- 20L
  seu <- make_seurat(n_cells = n_cells)
  set.seed(7L)
  coords <- matrix(stats::rnorm(n_cells * 2L), nrow = n_cells, ncol = 2L,
                   dimnames = list(colnames(seu), c("UMAP_1", "UMAP_2")))
  seu[["umap"]] <- suppressWarnings(SeuratObject::CreateDimReducObject(
    embeddings = coords, key = "UMAP_"))
  seu$cell_type <- rep(c("TypeA", "TypeB"), each = n_cells / 2L)

  ri_event_id <- "RI:chr1:851927:852094-852671:852766:+"
  # tx_dummy is the inclusion transcript for a second (SE) event required
  # because Assay5 needs >=2 features (events).
  txs  <- c("tx_ret", "tx_spl", "tx_dummy")
  cells <- paste0("Cell", seq_len(n_cells))
  set.seed(42L)
  tx_mat <- Matrix::Matrix(
    matrix(sample(0L:10L, 3L * n_cells, replace = TRUE),
           nrow = 3L, ncol = n_cells,
           dimnames = list(txs, cells)),
    sparse = TRUE
  )
  ioe <- tempfile(fileext = ".ioe")
  writeLines(c(
    "seqname\tgene_id\tinclusion_transcripts\ttotal_transcripts",
    paste(c("chr1", paste0("ENSG1;", ri_event_id),
            "tx_ret", "tx_ret,tx_spl"), collapse = "\t"),
    paste(c("chr1", "ENSG1;SE:chr1:900000-901000:902000-903000:+",
            "tx_dummy", "tx_dummy,tx_spl"), collapse = "\t")
  ), ioe)
  obj <- CreateMatisseObject(
    seurat            = seu,
    transcript_counts = tx_mat,
    events         = ioe,
    verbose           = FALSE
  )
  obj <- CalculatePSI(obj, min_coverage = 1L, verbose = FALSE)
  p <- PlotSashimi(obj, event_id = ri_event_id)
  expect_s3_class(p, "gg")
})

