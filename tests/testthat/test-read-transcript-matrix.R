# ---------------------------------------------------------------------------
# Tests for ReadTranscriptMatrix (long-read transcript triplet -> tx x cells)
# Reuses transcript fixtures from helper-fixtures.R.
# ---------------------------------------------------------------------------

# Write a (transcripts x cells) sparse matrix as an mtx triplet. `as_cells_x_tx`
# transposes on disk (Bagpiper-style). `gz` gzips the three files.
write_triplet <- function(m, dir = tempfile(), as_cells_x_tx = FALSE,
                          gz = FALSE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  disk <- if (as_cells_x_tx) Matrix::t(m) else m
  Matrix::writeMM(disk, file.path(dir, "matrix.mtx"))
  # features = transcript IDs, barcodes = cells (independent of disk layout)
  writeLines(rownames(m), file.path(dir, "features.tsv"))
  writeLines(colnames(m), file.path(dir, "barcodes.tsv"))
  if (gz) for (b in c("matrix.mtx", "barcodes.tsv", "features.tsv")) {
    f <- file.path(dir, b)
    con <- gzfile(paste0(f, ".gz"), "wb")
    writeBin(readBin(f, "raw", file.size(f)), con); close(con)
    file.remove(f)
  }
  dir
}

test_that("tx_x_cells triplet round-trips with dimnames (auto)", {
  m <- make_se_transcript_counts()                 # 5 tx x 20 cells
  d <- write_triplet(m)
  out <- ReadTranscriptMatrix(d, verbose = FALSE)
  expect_s4_class(out, "CsparseMatrix")
  expect_equal(rownames(out), rownames(m))
  expect_equal(colnames(out), colnames(m))
  expect_equal(as.matrix(out), as.matrix(m))
})

test_that("cells_x_tx source is auto-detected and transposed back", {
  m <- make_se_transcript_counts()
  d <- write_triplet(m, as_cells_x_tx = TRUE)       # Bagpiper layout
  out <- ReadTranscriptMatrix(d, verbose = FALSE)
  expect_equal(dim(out), dim(m))                    # transcripts x cells
  expect_equal(as.matrix(out), as.matrix(m))
})

test_that("gzipped triplet is read transparently", {
  m <- make_se_transcript_counts()
  d <- write_triplet(m, as_cells_x_tx = TRUE, gz = TRUE)
  out <- ReadTranscriptMatrix(d, verbose = FALSE)
  expect_equal(as.matrix(out), as.matrix(m))
})

test_that("explicit orientation mismatch errors clearly", {
  m <- make_se_transcript_counts()
  d <- write_triplet(m)                             # on-disk tx x cells
  expect_error(
    ReadTranscriptMatrix(d, orientation = "cells_x_tx", verbose = FALSE),
    "cells_x_tx")
})

test_that("strip_version collapses and sums versioned duplicates", {
  m <- make_se_transcript_counts(n_cells = 6L)
  rownames(m) <- c("ENST001.1", "ENST001.2", "ENST002.1",
                   "ENST003.4", "ENST004.2")
  d <- write_triplet(m)
  expect_warning(
    out <- ReadTranscriptMatrix(d, strip_version = TRUE, verbose = FALSE),
    "collapsed")
  expect_true("ENST001" %in% rownames(out))
  expect_equal(nrow(out), 4L)                       # 001.1 + 001.2 merged
  expect_equal(as.numeric(out["ENST001", ]),
               as.numeric(m["ENST001.1", ] + m["ENST001.2", ]))
})

test_that("strip_version = FALSE keeps version suffixes", {
  m <- make_se_transcript_counts(n_cells = 4L)
  rownames(m) <- c("ENST001.1", "ENST001.2", "ENST002.1",
                   "ENST003.4", "ENST004.2")
  d <- write_triplet(m)
  out <- ReadTranscriptMatrix(d, strip_version = FALSE, verbose = FALSE)
  expect_true(all(c("ENST001.1", "ENST001.2") %in% rownames(out)))
})

test_that("cells subsetting and zero-overlap behaviour", {
  m <- make_se_transcript_counts()
  d <- write_triplet(m)
  sub <- ReadTranscriptMatrix(d, cells = colnames(m)[1:5], verbose = FALSE)
  expect_equal(colnames(sub), colnames(m)[1:5])
  expect_error(ReadTranscriptMatrix(d, cells = "ZZZ", verbose = FALSE),
               "None of the requested")
})

test_that("missing dir / files error clearly", {
  expect_error(ReadTranscriptMatrix("/no/such/dir", verbose = FALSE),
               "does not exist")
  e <- tempfile(); dir.create(e)
  expect_error(ReadTranscriptMatrix(e, verbose = FALSE), "Could not find")
})

test_that("feeds CreateMatisseObject transcript mode end to end", {
  skip_if_not_installed("Seurat")
  m   <- make_se_transcript_counts(n_cells = 20L)
  d   <- write_triplet(m, as_cells_x_tx = TRUE)
  txc <- ReadTranscriptMatrix(d, verbose = FALSE)
  seu <- make_seurat(n_cells = 20L)                 # cells Cell1..Cell20
  obj <- CreateMatisseObject(seu, transcript_counts = txc,
                             events = make_se_ioe_file(),
                             min_coverage = 1L, verbose = FALSE)
  expect_s4_class(obj, "MatisseObject")
  expect_false(is.null(GetPSI(obj)))
})
