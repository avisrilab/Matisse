test_that("RunPCA dispatches on MatisseObject and returns MatisseObject", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  # NormalizeData / FindVariableFeatures / ScaleData are no longer Matisse
  # dispatchers; reach them via the embedded Seurat object instead.
  seu <- Seurat::NormalizeData(GetSeurat(obj), verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "dispersion",
                                       verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  obj@seurat <- seu
  result <- Seurat::RunPCA(obj, npcs = 5L, verbose = FALSE)
  expect_s4_class(result, "MatisseObject")
  expect_true("pca" %in% SeuratObject::Reductions(GetSeurat(result)))
})

test_that("AddMetaData dispatches on MatisseObject and returns MatisseObject", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  meta <- data.frame(batch = rep("B1", ncol(GetSeurat(obj))),
                     row.names = colnames(GetSeurat(obj)))
  result <- SeuratObject::AddMetaData(obj, metadata = meta)
  expect_s4_class(result, "MatisseObject")
  expect_true("batch" %in% colnames(result@seurat@meta.data))
})

test_that("FindMarkers dispatches on MatisseObject and returns a data frame", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  obj@seurat <- Seurat::NormalizeData(GetSeurat(obj), verbose = FALSE)
  # Set identities so both groups have cells
  SeuratObject::Idents(obj@seurat) <- rep(c("A", "B"), length.out = ncol(obj@seurat))
  result <- Seurat::FindMarkers(obj, ident.1 = "A", ident.2 = "B", verbose = FALSE)
  expect_s3_class(result, "data.frame")
})

test_that("SCTransform dispatches on a junction-mode MatisseObject (regression for @mode bug)", {
  # Regression: SCTransform.MatisseObject previously read object@mode == "event",
  # but the slot is input.mode and the value is "junction" / "transcript".
  # The function crashed on every call. Untested before this fix.
  skip_if_not_installed("Seurat")
  obj    <- make_matisse_object()  # junction mode
  result <- suppressWarnings(Seurat::SCTransform(obj, verbose = FALSE))
  expect_s4_class(result, "MatisseObject")
  expect_true("SCT" %in% SeuratObject::Assays(GetSeurat(result)))
})

test_that("DefaultAssay getter returns the embedded Seurat object's default assay", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  expect_equal(SeuratObject::DefaultAssay(obj),
               SeuratObject::DefaultAssay(GetSeurat(obj)))
})

test_that("DefaultAssay<- setter updates the embedded Seurat object", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  available <- SeuratObject::Assays(GetSeurat(obj))
  target    <- setdiff(available, SeuratObject::DefaultAssay(obj))[1]
  skip_if(is.na(target))
  SeuratObject::DefaultAssay(obj) <- target
  expect_equal(SeuratObject::DefaultAssay(obj), target)
  expect_equal(SeuratObject::DefaultAssay(GetSeurat(obj)), target)
})

test_that("FindMarkers exposes Seurat's named args (.onLoad formals copy)", {
  # Regression: before .onLoad was added the dispatcher had formals
  # `(object, ...)`, so RStudio tab-completion never offered ident.1 /
  # ident.2 / group.by / etc. Asserting these names appear in formals()
  # locks in the parameter-hint behaviour.
  expected <- c("ident.1", "ident.2", "group.by", "subset.ident",
                "assay", "latent.vars")
  fn <- get("FindMarkers.MatisseObject", envir = asNamespace("Matisse"))
  expect_true(all(expected %in% names(formals(fn))))
})

test_that("SCTransform dispatches on a transcript-mode MatisseObject", {
  skip_if_not_installed("Seurat")
  obj    <- make_matisse_from_transcripts()
  result <- suppressWarnings(Seurat::SCTransform(obj, verbose = FALSE))
  expect_s4_class(result, "MatisseObject")
  expect_true("SCT" %in% SeuratObject::Assays(GetSeurat(result)))
})
