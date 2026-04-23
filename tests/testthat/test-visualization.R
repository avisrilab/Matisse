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

test_that("PlotHeatmap: errors if PSI matrix is NULL", {
  obj <- make_matisse_object()
  expect_error(PlotHeatmap(obj), regexp = "PSI assay is NULL")
})

# ============================================================
# PlotCoverage
# ============================================================

test_that("PlotCoverage: returns a ggplot for a known gene", {
  obj <- make_matisse_with_umap()
  p   <- PlotCoverage(obj, gene = "gene1")
  expect_s3_class(p, "gg")
})

test_that("PlotCoverage: log_scale = TRUE still returns a ggplot", {
  obj <- make_matisse_with_umap()
  p   <- PlotCoverage(obj, gene = "gene1", log_scale = TRUE)
  expect_s3_class(p, "gg")
})

test_that("PlotCoverage: restricting to a cell subset returns a ggplot", {
  obj   <- make_matisse_with_umap()
  cells <- paste0("Cell", 1:5)
  p     <- PlotCoverage(obj, gene = "gene1", cells = cells)
  expect_s3_class(p, "gg")
})

test_that("PlotCoverage: bars are ordered by genomic start position", {
  obj <- make_matisse_with_umap()
  p   <- PlotCoverage(obj, gene = "gene1")
  # The x-axis levels should be in ascending start-position order
  jxn_levels  <- levels(p$data$junction)
  jd          <- GetJunctionData(obj)
  ordered_jxn <- jd$junction_id[order(jd$start)]
  ordered_jxn <- intersect(ordered_jxn, jxn_levels)   # keep only plotted junctions
  expect_equal(jxn_levels, ordered_jxn)
})

test_that("PlotCoverage: errors if no junction assay present", {
  skip_if_not_installed("Seurat")
  # Build an object without junction counts
  seu <- make_seurat()
  obj <- CreateMatisseObject(seurat = seu, verbose = FALSE)
  expect_error(PlotCoverage(obj, gene = "gene1"),
               regexp = "No junction assay found")
})

test_that("PlotCoverage: errors if junction_data is empty", {
  skip_if_not_installed("Seurat")
  seu <- make_seurat()
  jxn <- make_junction_counts()
  obj <- CreateMatisseObject(seurat = seu, junction_counts = jxn, verbose = FALSE)
  # junction_data is empty by default
  expect_error(PlotCoverage(obj, gene = "gene1"),
               regexp = "junction_data is empty")
})

test_that("PlotCoverage: errors for a gene not present in junction_data", {
  obj <- make_matisse_with_umap()
  expect_error(PlotCoverage(obj, gene = "no_such_gene"),
               regexp = "No junctions found")
})

# ============================================================
# CoveragePlot
# ============================================================

test_that("CoveragePlot: returns a ggplot for junction-mode object", {
  obj <- make_matisse_short_read()
  p   <- CoveragePlot(obj, event_id = "SE:chr1:1201-2999:3201-4999:+")
  expect_s3_class(p, "gg")
})

test_that("CoveragePlot: returns a ggplot for event-mode object", {
  obj <- make_matisse_long_read()
  p   <- CoveragePlot(obj, event_id = "SE:chr1:1201-2999:3201-4999:+")
  expect_s3_class(p, "gg")
})

test_that("CoveragePlot: group_by produces a faceted plot", {
  obj <- make_matisse_short_read()
  p   <- CoveragePlot(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                      group_by = "cell_type")
  expect_s3_class(p, "gg")
  expect_true(!is.null(p$facet))
})

test_that("CoveragePlot: arc_scale = 'linear' is accepted", {
  obj <- make_matisse_short_read()
  p   <- CoveragePlot(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                      arc_scale = "linear")
  expect_s3_class(p, "gg")
})

test_that("CoveragePlot: arc_scale = 'log' is accepted", {
  obj <- make_matisse_short_read()
  p   <- CoveragePlot(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                      arc_scale = "log")
  expect_s3_class(p, "gg")
})

test_that("CoveragePlot: custom title is applied", {
  obj <- make_matisse_short_read()
  p   <- CoveragePlot(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                      title = "My SE event")
  expect_equal(p$labels$title, "My SE event")
})

test_that("CoveragePlot: cell subset is accepted without error", {
  obj   <- make_matisse_short_read()
  cells <- paste0("Cell", 1:5)
  expect_no_error(
    CoveragePlot(obj, event_id = "SE:chr1:1201-2999:3201-4999:+",
                 cells = cells)
  )
})

test_that("CoveragePlot: errors for unknown event_id", {
  obj <- make_matisse_short_read()
  expect_error(
    CoveragePlot(obj, event_id = "SE:chr99:0-1:2-3:+"),
    regexp = "not found in event_data"
  )
})

test_that("CoveragePlot: errors in event mode for non-SE event_id format", {
  obj <- make_matisse_long_read()
  # Manually mangle the event_id to a non-SE type
  obj@event_data$event_id <- "A3:chr1:1201-2999:3201-4999:+"
  expect_error(
    CoveragePlot(obj, event_id = "A3:chr1:1201-2999:3201-4999:+"),
    regexp = "SE event"
  )
})

# ============================================================
# PlotQCMetrics
# ============================================================

test_that("PlotQCMetrics: returns a ggplot after ComputeIsoformQC", {
  obj <- make_matisse_with_qc()
  p   <- PlotQCMetrics(obj)
  expect_s3_class(p, "gg")
})

test_that("PlotQCMetrics: default features include all numeric QC columns", {
  obj      <- make_matisse_with_qc()
  meta     <- MatisseMeta(obj)
  num_cols <- colnames(meta)[vapply(meta, is.numeric, logical(1))]
  p        <- PlotQCMetrics(obj)
  # Facet labels should be a subset of numeric metadata columns
  facet_vals <- unique(p$data$metric)
  expect_true(all(facet_vals %in% num_cols))
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

test_that("PlotQCMetrics: errors for a feature not in metadata", {
  obj <- make_matisse_with_qc()
  expect_error(PlotQCMetrics(obj, features = "no_such_metric"),
               regexp = "not found")
})
