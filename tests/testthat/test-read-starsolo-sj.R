# ---------------------------------------------------------------------------
# Tests for ReadSTARsoloSJ (STARsolo SJ matrix -> cells x junctions)
# ---------------------------------------------------------------------------

# Build a minimal STARsolo SJ/raw directory.
#   features: 3 junctions, strand codes 1 (+), 2 (-), 0 (undefined)
#   matrix:   junctions x cells
make_sj_dir <- function(dir = tempfile()) {
  dir.create(dir, recursive = TRUE)
  # col4 = STAR strand code, col5 = intron motif code
  feat <- rbind(
    c("chr1", "100", "200", "1", "1", "1", "5", "0", "30"),  # + (code)
    c("chr1", "300", "400", "2", "2", "1", "3", "0", "28"),  # - (code)
    c("chr2", "500", "600", "0", "0", "0", "1", "0", "20"),  # undefined
    c("chr3", "700", "800", "0", "3", "0", "2", "0", "22")   # undef + GC/AG
  )
  utils::write.table(feat, file.path(dir, "features.tsv"),
                     sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  writeLines(c("CellA", "CellB", "CellC"),
             file.path(dir, "barcodes.tsv"))
  m <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3, 4), j = c(1, 2, 2, 3, 1), x = c(4, 7, 2, 9, 6),
    dims = c(4L, 3L))
  Matrix::writeMM(m, file.path(dir, "matrix.mtx"))
  dir
}

test_that("returns cells x junctions with chr-start-end-strand IDs", {
  d   <- make_sj_dir()
  jxn <- ReadSTARsoloSJ(d, verbose = FALSE)
  expect_s4_class(jxn, "CsparseMatrix")
  expect_equal(rownames(jxn), c("CellA", "CellB", "CellC"))
  # junction 3 is strand-undefined with non-canonical motif -> "*";
  # junction 4 is strand-undefined but GC/AG (motif 3) -> inferred "+"
  expect_setequal(colnames(jxn),
                  c("chr1-100-200-+", "chr1-300-400--",
                    "chr2-500-600-*", "chr3-700-800-+"))
  expect_equal(jxn["CellB", "chr1-100-200-+"], 7)
  expect_equal(jxn["CellB", "chr1-300-400--"], 2)
})

test_that("infer_strand recovers strand from intron motif", {
  d  <- make_sj_dir()
  on <- ReadSTARsoloSJ(d, verbose = FALSE)
  expect_true("chr3-700-800-+" %in% colnames(on))   # motif 3 -> +
  off <- ReadSTARsoloSJ(d, infer_strand = FALSE, verbose = FALSE)
  expect_true("chr3-700-800-*" %in% colnames(off))  # left undefined
})

test_that("strand_map controls strand symbols", {
  d   <- make_sj_dir()
  jxn <- ReadSTARsoloSJ(d, strand_map = c("0" = ".", "1" = "+", "2" = "-"),
                        infer_strand = FALSE, verbose = FALSE)
  expect_true("chr2-500-600-." %in% colnames(jxn))
})

test_that("cells subsetting keeps only requested barcodes", {
  d   <- make_sj_dir()
  jxn <- ReadSTARsoloSJ(d, cells = c("CellA", "CellC"), verbose = FALSE)
  expect_equal(rownames(jxn), c("CellA", "CellC"))
})

test_that("no overlapping cells errors clearly", {
  d <- make_sj_dir()
  expect_error(ReadSTARsoloSJ(d, cells = c("ZZZ"), verbose = FALSE),
               "None of the requested")
})

test_that("missing directory / files error clearly", {
  expect_error(ReadSTARsoloSJ("/no/such/dir", verbose = FALSE),
               "does not exist")
  d <- tempfile(); dir.create(d)
  expect_error(ReadSTARsoloSJ(d, verbose = FALSE), "Could not find")
})

test_that("counts feed CreateMatisseObject junction mode end to end", {
  skip_if_not_installed("Seurat")
  d   <- make_sj_dir()
  jxn <- ReadSTARsoloSJ(d, verbose = FALSE)
  # >=2 events (Assay5 requires >=2 features) using real SJ junction IDs
  ev <- data.frame(
    event_id           = c("EV1", "EV2"),
    gene_id            = c("G1", "G2"),
    chr                = c("chr1", "chr1"),
    strand             = c("+", "-"),
    event_type         = c("SE", "SE"),
    inclusion_features = c("chr1-100-200-+", "chr1-300-400--"),
    exclusion_features = c("chr1-300-400--", "chr2-500-600-*"),
    stringsAsFactors   = FALSE)
  gene <- matrix(rpois(60, 5), nrow = 20,
                 dimnames = list(paste0("G", 1:20),
                                 c("CellA", "CellB", "CellC")))
  seu  <- Seurat::CreateSeuratObject(gene)
  obj  <- CreateMatisseObject(seu, junction_counts = jxn, events = ev,
                              min_coverage = 1L, verbose = FALSE)
  expect_s4_class(obj, "MatisseObject")
})
