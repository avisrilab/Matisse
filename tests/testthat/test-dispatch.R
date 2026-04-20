test_that("NormalizeData dispatches on MatisseObject and returns MatisseObject", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  result <- Seurat::NormalizeData(obj, verbose = FALSE)
  expect_s4_class(result, "MatisseObject")
})

test_that("FindVariableFeatures dispatches on MatisseObject and returns MatisseObject", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  result <- Seurat::FindVariableFeatures(obj, selection.method = "dispersion",
                                         verbose = FALSE)
  expect_s4_class(result, "MatisseObject")
})

test_that("ScaleData dispatches on MatisseObject and returns MatisseObject", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "dispersion",
                                      verbose = FALSE)
  result <- Seurat::ScaleData(obj, verbose = FALSE)
  expect_s4_class(result, "MatisseObject")
})

test_that("RunPCA dispatches on MatisseObject and returns MatisseObject", {
  skip_if_not_installed("Seurat")
  obj <- make_matisse_object()
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "dispersion",
                                      verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
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
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  # Set identities so both groups have cells
  SeuratObject::Idents(obj@seurat) <- rep(c("A", "B"), length.out = ncol(obj@seurat))
  result <- Seurat::FindMarkers(obj, ident.1 = "A", ident.2 = "B", verbose = FALSE)
  expect_s3_class(result, "data.frame")
})
