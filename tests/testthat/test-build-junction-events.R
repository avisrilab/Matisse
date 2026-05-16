# ---------------------------------------------------------------------------
# Tests for BuildJunctionEvents (SUPPA2 event_id -> junction-ID adapter)
# ---------------------------------------------------------------------------

write_ioe <- function(rows, file = tempfile(fileext = ".ioe")) {
  writeLines(c(
    "seqname\tgene_id\tinclusion_transcripts\ttotal_transcripts",
    rows), file)
  file
}

test_that("SE: inclusion = two flanking junctions, exclusion = skipping", {
  f <- write_ioe(paste(c("chr1",
    "ENSG1;SE:chr1:100-200:300-400:+", "txA", "txA,txB"), collapse = "\t"))
  ev <- BuildJunctionEvents(f, junction_universe = NULL, verbose = FALSE)
  expect_equal(nrow(ev), 1L)
  expect_equal(ev$event_type, "SE")
  expect_equal(ev$inclusion_features, "chr1-101-199-+;chr1-301-399-+")
  expect_equal(ev$exclusion_features, "chr1-101-399-+")
})

test_that("A3/A5: one discriminating junction each side", {
  f <- write_ioe(c(
    paste(c("chr1", "ENSG1;A3:chr1:50-100:50-120:+",
            "txA", "txA,txB"), collapse = "\t"),
    paste(c("chr1", "ENSG2;A5:chr1:200-300:220-300:-",
            "txC", "txC,txD"), collapse = "\t")))
  ev <- BuildJunctionEvents(f, verbose = FALSE)
  a3 <- ev[ev$event_type == "A3", ]
  expect_equal(a3$inclusion_features, "chr1-51-99-+")
  expect_equal(a3$exclusion_features, "chr1-51-119-+")
  a5 <- ev[ev$event_type == "A5", ]
  expect_equal(a5$inclusion_features, "chr1-201-299--")
  expect_equal(a5$exclusion_features, "chr1-221-299--")
})

test_that("MX: two junctions per mutually-exclusive exon", {
  f <- write_ioe(paste(c("chr1",
    "ENSG1;MX:chr1:10-20:30-40:10-50:60-70:+", "txA", "txA,txB"),
    collapse = "\t"))
  ev <- BuildJunctionEvents(f, verbose = FALSE)
  expect_equal(ev$inclusion_features, "chr1-11-19-+;chr1-31-39-+")
  expect_equal(ev$exclusion_features, "chr1-11-49-+;chr1-61-69-+")
})

test_that("AF/AL: pair blocks become the discriminating junctions", {
  f <- write_ioe(c(
    paste(c("chr1", "ENSG1;AF:chr1:5:10-100:55:60-100:+",
            "txA", "txA,txB"), collapse = "\t"),
    paste(c("chr1", "ENSG2;AL:chr1:10-50:60:10-80:90:+",
            "txC", "txC,txD"), collapse = "\t")))
  ev <- BuildJunctionEvents(f, verbose = FALSE)
  af <- ev[ev$event_type == "AF", ]
  expect_equal(af$inclusion_features, "chr1-11-99-+")
  expect_equal(af$exclusion_features, "chr1-61-99-+")
  al <- ev[ev$event_type == "AL", ]
  expect_equal(al$inclusion_features, "chr1-11-49-+")
  expect_equal(al$exclusion_features, "chr1-11-79-+")
})

test_that("RI is dropped with a warning (not junction-quantifiable)", {
  f <- write_ioe(c(
    paste(c("chr1", "ENSG1;SE:chr1:100-200:300-400:+",
            "txA", "txA,txB"), collapse = "\t"),
    paste(c("chr1", "ENSG2;RI:chr1:100:200-300:400:+",
            "txC", "txC,txD"), collapse = "\t")))
  expect_warning(
    ev <- BuildJunctionEvents(f, event_types = c("SE", "RI"),
                              verbose = FALSE),
    "junction-quantifiable")
  expect_equal(nrow(ev), 1L)
  expect_equal(ev$event_type, "SE")
})

test_that("explicit integer intron_offset is honoured", {
  f <- write_ioe(paste(c("chr1",
    "ENSG1;SE:chr1:100-200:300-400:+", "txA", "txA,txB"), collapse = "\t"))
  ev <- BuildJunctionEvents(f, intron_offset = c(0L, 0L), verbose = FALSE)
  expect_equal(ev$inclusion_features, "chr1-100-200-+;chr1-300-400-+")
  expect_equal(ev$exclusion_features, "chr1-100-400-+")
})

test_that("auto intron_offset calibrates to the universe", {
  f <- write_ioe(paste(c("chr1",
    "ENSG1;SE:chr1:100-200:300-400:+", "txA", "txA,txB"), collapse = "\t"))
  # Universe encodes the c(0,0) convention; auto must pick it.
  uni <- c("chr1-100-200-+", "chr1-300-400-+", "chr1-100-400-+")
  ev  <- BuildJunctionEvents(f, junction_universe = uni,
                             intron_offset = "auto", verbose = FALSE)
  expect_equal(ev$inclusion_features, "chr1-100-200-+;chr1-300-400-+")
})

test_that("missing IOE files error clearly", {
  expect_error(BuildJunctionEvents("does-not-exist.ioe", verbose = FALSE),
               "not found")
})
