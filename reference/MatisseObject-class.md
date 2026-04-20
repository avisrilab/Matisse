# The MatisseObject S4 class

The central data structure for Matisse. It wraps a
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
object and augments it with isoform-resolved layers. Two fixed assays
live inside the embedded Seurat object:

- `"transcript"`:

  A standard `Assay5` holding raw transcript-level counts (transcripts ×
  cells). Created by
  [`CreateMatisseObjectFromTranscripts`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObjectFromTranscripts.md).

- `"psi"`:

  An `Assay5` holding PSI values in the `"data"` layer, inclusion counts
  in `"counts"`, and exclusion counts in `"exclusion"` (all features ×
  cells). Created by
  [`CalculatePSI`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
  or
  [`CreateMatisseObjectFromTranscripts`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObjectFromTranscripts.md).

Raw junction counts and splice-event annotations are kept as slots on
the MatisseObject itself.

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

## Functions

- `show(MatisseObject)`: Display a summary of a `MatisseObject`.

- `x[[i`: Access isoform metadata columns or Seurat slots via `[[`.
  Checks `isoform_metadata` first; falls back to the embedded Seurat
  object.

## Slots

- `seurat`:

  A `Seurat` object carrying gene-level expression, dimensionality
  reductions, cell metadata, and the `"transcript"` / `"psi"` assays.

- `junction_counts`:

  A sparse matrix (dgCMatrix, cells × junctions) of raw per-junction
  read counts.

- `event_data`:

  A `data.frame` with one row per splice event. Required columns:
  `event_id`, `gene_id`, `chr`, `strand`, `event_type`,
  `inclusion_junctions`, `exclusion_junctions`.

- `junction_data`:

  A `data.frame` with one row per junction. Required columns:
  `junction_id`, `chr`, `start`, `end`, `strand`, `gene_id`.

- `isoform_metadata`:

  A `data.frame` of per-cell isoform QC metrics. Rownames correspond to
  cell barcodes.

- `version`:

  Character string recording the Matisse version used to create the
  object.

- `misc`:

  Named list for user-defined extra data.
