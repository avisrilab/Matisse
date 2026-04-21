# The MatisseObject S4 class

The central data structure for Matisse. It wraps a
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
object and adds isoform-resolved splicing layers. All per-cell data —
junction counts, PSI values, transcript counts, and QC metrics — live
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

Two operating modes are supported, set automatically at construction:

- `"junction"`:

  Short-read mode. Raw junction counts are stored as
  `Assay5("junction")` (junctions × cells). PSI is computed later by
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  and stored as `Assay5("psi")`.

- `"event"`:

  Long-read mode. Transcript counts (e.g. from Bagpiper or FLAMES) are
  stored as `Assay5("transcript")`. PSI is computed at construction time
  from SUPPA2 `.ioe` event definitions and stored as `Assay5("psi")`.

## Functions

- `show(MatisseObject)`: Display a summary of a `MatisseObject`.

- `x[[i`: Access cell metadata or Seurat slots via `[[`. Checks
  `seurat@meta.data` first; falls back to the embedded Seurat object
  (assays, reductions, etc.).

## Slots

- `seurat`:

  A `Seurat` object. Contains all per-cell data: gene expression, splice
  assays (`"junction"`, `"transcript"`, `"psi"`), cell metadata (QC
  metrics, cluster labels), and dimensionality reductions.

- `event_data`:

  A `data.frame` with one row per splice event. Required columns:
  `event_id`, `gene_id`, `chr`, `strand`, `event_type`,
  `inclusion_junctions`, `exclusion_junctions`.

- `junction_data`:

  A `data.frame` with one row per junction. Required columns:
  `junction_id`, `chr`, `start`, `end`, `strand`, `gene_id`.

- `mode`:

  Character. `"junction"` for short-read objects; `"event"` for
  long-read objects.

- `version`:

  Character string recording the Matisse version used to create the
  object.

- `misc`:

  Named list for user-defined extra data.
