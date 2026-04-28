# The MatisseObject S4 class

The central data structure for Matisse. It wraps a
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
object and adds isoform-resolved splicing layers. All per-cell data –
junction counts, PSI values, transcript counts, and QC metrics – live
inside the embedded Seurat object as named assays (`Assay5`) or cell
metadata (`meta.data`). Nothing is duplicated outside the Seurat object.

## Usage

``` r
# S4 method for class 'MatisseObject'
show(object)

# S4 method for class 'MatisseObject'
dim(x)

# S4 method for class 'MatisseObject,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MatisseObject,ANY,ANY'
x[[i, j, ...]]

# S4 method for class 'MatisseObject'
x$name
```

## Arguments

- object:

  A `MatisseObject`.

- x:

  A `MatisseObject`.

- i:

  Cell barcodes (character) or integer indices.

- j:

  Event IDs (character) or integer indices.

- ...:

  Ignored.

- drop:

  Ignored.

- name:

  Metadata column name or Seurat/Signac function name.

## Details

Two input modes are supported, set automatically at construction:

- `"junction"`:

  Short-read mode. Raw junction counts (from STARsolo) are stored as
  `Assay5("isoform")` (junctions x cells). Call
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  to compute PSI values.

- `"transcript"`:

  Long-read mode. Transcript isoform counts (from Bagpiper, FLAMES, or
  LIQA) are stored as `Assay5("isoform")` (transcripts x cells). Call
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  to compute PSI values.

After
[`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md),
PSI values are stored as `Assay5("psi")` (splice events x cells).

## Functions

- `show(MatisseObject)`: Display a summary of a `MatisseObject`.

- `x[[i`: Access cell metadata or Seurat slots via `[[`. Checks
  `seurat@meta.data` first; falls back to the embedded Seurat object
  (assays, reductions, etc.).

## Slots

- `seurat`:

  A `Seurat` object. Contains all per-cell data: gene expression, splice
  assays (`"isoform"`, `"psi"`), cell metadata (QC metrics, cluster
  labels), and dimensionality reductions.

- `input.mode`:

  Character. `"junction"` for short-read (STARsolo junction counts)
  objects; `"transcript"` for long-read (Bagpiper / FLAMES / LIQA
  transcript counts) objects.

- `version`:

  Character string recording the Matisse version used to create the
  object.

- `misc`:

  Named list for user-defined extra data.
