# Run PCA on a MatisseObject

Runs [`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html) on
the embedded Seurat object and returns the updated `MatisseObject`. The
PCA result is stored inside the Seurat object and accessible via
`GetSeurat(obj)`.

## Usage

``` r
# S3 method for class 'MatisseObject'
RunPCA(
  object,
  assay = NULL,
  features = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42,
  ...
)
```

## Arguments

- object:

  A `MatisseObject`.

- assay:

  Name of Assay PCA is being run on

- features:

  Features to compute PCA on. If features=NULL, PCA will be run using
  the variable features for the Assay. Note that the features must be
  present in the scaled data. Any requested features that are not scaled
  or have 0 variance will be dropped, and the PCA will be run using the
  remaining features.

- npcs:

  Total Number of PCs to compute and store (50 by default)

- rev.pca:

  By default computes the PCA on the cell x gene matrix. Setting to true
  will compute it on gene x cell matrix.

- weight.by.var:

  Weight the cell embeddings by the variance of each PC (weights the
  gene loadings if rev.pca is TRUE)

- verbose:

  Print the top genes associated with high/low loadings for the PCs

- ndims.print:

  PCs to print genes for

- nfeatures.print:

  Number of genes to print for each PC

- reduction.name:

  dimensional reduction name, pca by default

- reduction.key:

  dimensional reduction key, specifies the string before the number for
  the dimension names. PC by default

- seed.use:

  Set a random seed. By default, sets the seed to 42. Setting NULL will
  not set a seed.

- ...:

  Additional arguments forwarded to
  [`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html) (e.g.
  `assay`, `npcs`, `features`).

## Value

The updated `MatisseObject` with a `"pca"` reduction.

## Details

Typical usage after
[`SCTransform.MatisseObject`](https://avisrilab.org/Matisse/reference/SCTransform.MatisseObject.md):


    obj <- RunPCA(obj, assay = "SCT", npcs = 50)

## See also

[`RunUMAP.MatisseObject`](https://avisrilab.org/Matisse/reference/RunUMAP.MatisseObject.md),
[`SCTransform.MatisseObject`](https://avisrilab.org/Matisse/reference/SCTransform.MatisseObject.md)
