# ---------------------------------------------------------------------------
# Tests for BuildSimpleEvents
# ---------------------------------------------------------------------------

test_that("BuildSimpleEvents: returns a data.frame with all required columns", {
  junctions <- make_junction_data()
  result    <- BuildSimpleEvents(junctions)
  expect_s3_class(result, "data.frame")
  required <- c("event_id", "gene_id", "chr", "strand",
                "event_type", "inclusion_junctions", "exclusion_junctions")
  expect_true(all(required %in% colnames(result)))
})

test_that("BuildSimpleEvents: produces exactly one row per junction", {
  junctions <- make_junction_data()   # 6 junctions
  result    <- BuildSimpleEvents(junctions)
  expect_equal(nrow(result), nrow(junctions))
})

test_that("BuildSimpleEvents: event_id is formatted as gene_id:junction_id", {
  junctions <- make_junction_data()
  result    <- BuildSimpleEvents(junctions)
  # Every event_id must start with the gene_id and contain the junction_id
  expect_true(all(startsWith(result$event_id, "gene1:")))
  jxn_parts <- sub("^gene1:", "", result$event_id)
  expect_setequal(jxn_parts, junctions$junction_id)
})

test_that("BuildSimpleEvents: inclusion_junctions is the focal junction only", {
  junctions <- make_junction_data()
  result    <- BuildSimpleEvents(junctions)
  # No semicolons — each inclusion field is a single junction
  expect_false(any(grepl(";", result$inclusion_junctions)))
  expect_true(all(result$inclusion_junctions %in% junctions$junction_id))
})

test_that("BuildSimpleEvents: exclusion_junctions contains every other junction for the gene", {
  junctions <- make_junction_data()   # 6 junctions, all gene1
  result    <- BuildSimpleEvents(junctions)
  for (i in seq_len(nrow(result))) {
    focal   <- result$inclusion_junctions[i]
    exc_str <- result$exclusion_junctions[i]
    if (nchar(exc_str) == 0) next
    exc_ids <- strsplit(exc_str, ";", fixed = TRUE)[[1]]
    expect_false(focal %in% exc_ids,
                 info = paste("focal junction", focal, "appeared in its own exclusion set"))
    expect_setequal(
      sort(c(focal, exc_ids)),
      sort(junctions$junction_id),
      info  = paste("event for", focal, "- union of inc+exc differs from all junctions")
    )
  }
})

test_that("BuildSimpleEvents: single-junction gene produces empty exclusion_junctions", {
  single <- data.frame(
    junction_id      = "jxnOnly",
    gene_id          = "geneX",
    chr              = "chr2",
    strand           = "-",
    stringsAsFactors = FALSE
  )
  result <- BuildSimpleEvents(single)
  expect_equal(nrow(result), 1L)
  expect_equal(result$exclusion_junctions, "")
})

test_that("BuildSimpleEvents: junctions from different genes do not cross-contaminate", {
  two_genes <- data.frame(
    junction_id      = c("jA", "jB", "jC", "jD"),
    gene_id          = c("g1",  "g1",  "g2",  "g2"),
    chr              = rep("chr1", 4L),
    strand           = rep("+",   4L),
    stringsAsFactors = FALSE
  )
  result <- BuildSimpleEvents(two_genes)
  expect_equal(nrow(result), 4L)

  # jC's exclusion should be jD only — not jA or jB
  row_jC  <- result[result$inclusion_junctions == "jC", ]
  exc_jC  <- strsplit(row_jC$exclusion_junctions, ";", fixed = TRUE)[[1]]
  expect_false("jA" %in% exc_jC)
  expect_false("jB" %in% exc_jC)
  expect_true("jD"  %in% exc_jC)
})

test_that("BuildSimpleEvents: chr and strand are passed through from the junction table", {
  junctions <- make_junction_data()
  result    <- BuildSimpleEvents(junctions)
  expect_true(all(result$chr    == "chr1"))
  expect_true(all(result$strand == "+"))
})

test_that("BuildSimpleEvents: event_type column is always 'simple'", {
  junctions <- make_junction_data()
  result    <- BuildSimpleEvents(junctions)
  expect_true(all(result$event_type == "simple"))
})

test_that("BuildSimpleEvents: gene_id column is preserved", {
  junctions <- make_junction_data()
  result    <- BuildSimpleEvents(junctions)
  expect_true(all(result$gene_id == "gene1"))
})

test_that("BuildSimpleEvents: result is usable as event_data in CreateMatisseObject", {
  skip_if_not_installed("Seurat")
  seu    <- make_seurat()
  jxn    <- make_junction_counts()
  events <- BuildSimpleEvents(make_junction_data())
  expect_no_error(
    CreateMatisseObject(
      seurat          = seu,
      junction_counts = jxn,
      event_data      = events,
      verbose         = FALSE
    )
  )
})

test_that("BuildSimpleEvents: errors for unsupported event_type", {
  junctions <- make_junction_data()
  expect_error(
    BuildSimpleEvents(junctions, event_type = "SE"),
    regexp = "Only event_type"
  )
})

test_that("BuildSimpleEvents: errors when required column 'junction_id' is missing", {
  bad <- data.frame(gene_id = "g1", stringsAsFactors = FALSE)
  expect_error(BuildSimpleEvents(bad), regexp = "junction_id")
})

test_that("BuildSimpleEvents: errors when required column 'gene_id' is missing", {
  bad <- data.frame(junction_id = "jxn1", stringsAsFactors = FALSE)
  expect_error(BuildSimpleEvents(bad), regexp = "gene_id")
})
