# ---------------------------------------------------------------------------
# Tests for visualization functions
# All Plot* functions return a ggplot object; we do not snapshot pixel output.
# ---------------------------------------------------------------------------

# ============================================================
# PlotPSIUMAP
# ============================================================

test_that("PlotPSIUMAP: returns a ggplot object", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIUMAP(obj, event_id = "SE-gene1-e2")
  expect_s3_class(p, "gg")
})

test_that("PlotPSIUMAP: custom title is applied", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIUMAP(obj, event_id = "SE-gene1-e2", title = "My Title")
  expect_equal(p$labels$title, "My Title")
})

test_that("PlotPSIUMAP: default title is the event_id", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIUMAP(obj, event_id = "SE-gene1-e2")
  expect_equal(p$labels$title, "SE-gene1-e2")
})

test_that("PlotPSIUMAP: data has the same number of rows as the cell count", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIUMAP(obj, event_id = "SE-gene1-e2")
  expect_equal(nrow(p$data), .n_cells(obj))
})

test_that("PlotPSIUMAP: PSI values are in [0, 1] or NA in the plot data", {
  obj  <- make_matisse_with_umap()
  p    <- PlotPSIUMAP(obj, event_id = "SE-gene1-e2")
  vals <- p$data$psi
  finite_vals <- vals[!is.na(vals)]
  expect_true(all(finite_vals >= 0 & finite_vals <= 1))
})

test_that("PlotPSIUMAP: errors if PSI matrix is NULL", {
  obj <- make_matisse_object()   # PSI not yet calculated
  expect_error(PlotPSIUMAP(obj, event_id = "SE-gene1-e2"),
               regexp = "PSI assay is NULL")
})

test_that("PlotPSIUMAP: errors for an event_id not present in the PSI matrix", {
  obj <- make_matisse_with_umap()
  expect_error(PlotPSIUMAP(obj, event_id = "nonexistent_event_XYZ"),
               regexp = "not found")
})

test_that("PlotPSIUMAP: errors if the requested reduction is absent from the Seurat object", {
  # make_matisse_object uses a seurat WITHOUT a umap reduction
  obj <- CalculatePSI(make_matisse_object(), min_coverage = 1L, verbose = FALSE)
  expect_error(PlotPSIUMAP(obj, event_id = "SE-gene1-e2"))
})

# ============================================================
# PlotPSIViolin
# ============================================================

test_that("PlotPSIViolin: returns a ggplot object", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIViolin(obj, event_id = "SE-gene1-e2", group_by = "cell_type")
  expect_s3_class(p, "gg")
})

test_that("PlotPSIViolin: default group_by 'seurat_clusters' is used when not specified", {
  skip_if_not_installed("Seurat")
  # Seurat always creates seurat_clusters (set to 0 by default) in CreateSeuratObject
  obj <- make_matisse_with_umap()
  expect_no_error(PlotPSIViolin(obj, event_id = "SE-gene1-e2"))
})

test_that("PlotPSIViolin: add_points = TRUE still returns a ggplot", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIViolin(obj, event_id = "SE-gene1-e2",
                        group_by = "cell_type", add_points = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotPSIViolin: custom colour vector is accepted", {
  obj    <- make_matisse_with_umap()
  colors <- c(TypeA = "#E41A1C", TypeB = "#377EB8")
  p      <- PlotPSIViolin(obj, event_id = "SE-gene1-e2",
                           group_by = "cell_type", colours = colors)
  expect_s3_class(p, "gg")
})

test_that("PlotPSIViolin: custom title is applied", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIViolin(obj, event_id = "SE-gene1-e2",
                        group_by = "cell_type", title = "Splicing switch")
  expect_equal(p$labels$title, "Splicing switch")
})

test_that("PlotPSIViolin: errors if PSI matrix is NULL", {
  obj <- make_matisse_object()
  expect_error(PlotPSIViolin(obj, event_id = "SE-gene1-e2"),
               regexp = "PSI assay is NULL")
})

test_that("PlotPSIViolin: errors for an unknown event_id", {
  obj <- make_matisse_with_umap()
  expect_error(PlotPSIViolin(obj, event_id = "ghost_event"),
               regexp = "not found")
})

test_that("PlotPSIViolin: errors for an unknown group_by column", {
  obj <- make_matisse_with_umap()
  expect_error(
    PlotPSIViolin(obj, event_id = "SE-gene1-e2", group_by = "no_such_col"),
    regexp = "not found"
  )
})

# ============================================================
# PlotPSIHeatmap
# ============================================================

test_that("PlotPSIHeatmap: returns a ggplot object", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIHeatmap(obj)
  expect_s3_class(p, "gg")
})

test_that("PlotPSIHeatmap: max_cells downsampling is respected", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIHeatmap(obj, max_cells = 5L)
  expect_lte(length(unique(p$data$cell)), 5L)
})

test_that("PlotPSIHeatmap: subsetting to specific events works", {
  obj <- make_matisse_with_umap()
  p   <- PlotPSIHeatmap(obj, events = "SE-gene1-e2")
  expect_equal(length(unique(p$data$event)), 1L)
})

test_that("PlotPSIHeatmap: subsetting to specific cells works", {
  obj   <- make_matisse_with_umap()
  cells <- paste0("Cell", 1:5)
  p     <- PlotPSIHeatmap(obj, cells = cells)
  expect_lte(length(unique(p$data$cell)), 5L)
})

test_that("PlotPSIHeatmap: group_by orders cells without error", {
  obj <- make_matisse_with_umap()
  expect_no_error(PlotPSIHeatmap(obj, group_by = "cell_type"))
})

test_that("PlotPSIHeatmap: warns on unknown event IDs (does not error)", {
  obj <- make_matisse_with_umap()
  expect_warning(
    PlotPSIHeatmap(obj, events = c("SE-gene1-e2", "bad_event_99")),
    regexp = "not found"
  )
})

test_that("PlotPSIHeatmap: errors if PSI matrix is NULL", {
  obj <- make_matisse_object()
  expect_error(PlotPSIHeatmap(obj), regexp = "PSI assay is NULL")
})

# ============================================================
# PlotJunctionCoverage
# ============================================================

test_that("PlotJunctionCoverage: returns a ggplot for a known gene", {
  obj <- make_matisse_with_umap()
  p   <- PlotJunctionCoverage(obj, gene = "gene1")
  expect_s3_class(p, "gg")
})

test_that("PlotJunctionCoverage: log_scale = TRUE still returns a ggplot", {
  obj <- make_matisse_with_umap()
  p   <- PlotJunctionCoverage(obj, gene = "gene1", log_scale = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotJunctionCoverage: restricting to a cell subset returns a ggplot", {
  obj   <- make_matisse_with_umap()
  cells <- paste0("Cell", 1:5)
  p     <- PlotJunctionCoverage(obj, gene = "gene1", cells = cells)
  expect_s3_class(p, "gg")
})

test_that("PlotJunctionCoverage: bars are ordered by genomic start position", {
  obj <- make_matisse_with_umap()
  p   <- PlotJunctionCoverage(obj, gene = "gene1")
  # The x-axis levels should be in ascending start-position order
  jxn_levels  <- levels(p$data$junction)
  jd          <- GetJunctionData(obj)
  ordered_jxn <- jd$junction_id[order(jd$start)]
  ordered_jxn <- intersect(ordered_jxn, jxn_levels)   # keep only plotted junctions
  expect_equal(jxn_levels, ordered_jxn)
})

test_that("PlotJunctionCoverage: errors if junction_counts slot is NULL", {
  skip_if_not_installed("Seurat")
  # Build an object without junction counts
  seu <- make_seurat()
  obj <- CreateMatisseObject(seurat = seu, verbose = FALSE)
  expect_error(PlotJunctionCoverage(obj, gene = "gene1"),
               regexp = "junction_counts slot is NULL")
})

test_that("PlotJunctionCoverage: errors if junction_data is empty", {
  skip_if_not_installed("Seurat")
  seu <- make_seurat()
  jxn <- make_junction_counts()
  obj <- CreateMatisseObject(seurat = seu, junction_counts = jxn, verbose = FALSE)
  # junction_data is empty by default
  expect_error(PlotJunctionCoverage(obj, gene = "gene1"),
               regexp = "junction_data is empty")
})

test_that("PlotJunctionCoverage: errors for a gene not present in junction_data", {
  obj <- make_matisse_with_umap()
  expect_error(PlotJunctionCoverage(obj, gene = "no_such_gene"),
               regexp = "No junctions found")
})

# ============================================================
# PlotQCMetrics
# ============================================================

test_that("PlotQCMetrics: returns a ggplot after ComputeIsoformQC", {
  obj <- make_matisse_with_qc()
  p   <- PlotQCMetrics(obj)
  expect_s3_class(p, "gg")
})

test_that("PlotQCMetrics: default features are all numeric columns", {
  obj      <- make_matisse_with_qc()
  meta     <- MatisseMeta(obj)
  num_cols <- colnames(meta)[vapply(meta, is.numeric, logical(1))]
  p        <- PlotQCMetrics(obj)
  # Facet labels should cover all numeric QC columns
  facet_vals <- unique(p$data$metric)
  expect_setequal(facet_vals, num_cols)
})

test_that("PlotQCMetrics: subset of features is accepted", {
  obj <- make_matisse_with_qc()
  p   <- PlotQCMetrics(obj, features = "mean_psi")
  expect_equal(unique(p$data$metric), "mean_psi")
})

test_that("PlotQCMetrics: group_by splits cells into multiple violin groups", {
  obj <- make_matisse_with_qc()
  p   <- PlotQCMetrics(obj, features = "mean_psi", group_by = "cell_type")
  expect_true(length(unique(p$data$group)) > 1L)
})

test_that("PlotQCMetrics: ncol is forwarded to the facet without error", {
  obj <- make_matisse_with_qc()
  p   <- PlotQCMetrics(obj, ncol = 1L)
  expect_s3_class(p, "gg")
})

test_that("PlotQCMetrics: errors for a feature not in isoform_metadata", {
  obj <- make_matisse_with_qc()
  expect_error(PlotQCMetrics(obj, features = "no_such_metric"),
               regexp = "not found")
})
