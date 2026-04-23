# Calculate PSI matrix from junction counts

Computes a Percent Spliced In (PSI) matrix for all splice events defined
in `event_data`. Only applies to objects in **junction mode**; in event
mode PSI is computed at construction time.

## Usage

``` r
CalculatePSI(
  object,
  events = NULL,
  min_coverage = 5L,
  na_fill = NA_real_,
  verbose = TRUE,
  ...
)

# S4 method for class 'MatisseObject'
CalculatePSI(
  object,
  events = NULL,
  min_coverage = 5L,
  na_fill = NA_real_,
  verbose = TRUE,
  ...
)

# S4 method for class 'ANY'
CalculatePSI(
  object,
  events = NULL,
  min_coverage = 5L,
  na_fill = NA_real_,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A
  [`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  in junction mode, or a sparse matrix (cells × junctions).

- events:

  When `object` is a matrix: a `data.frame` with columns `event_id`,
  `inclusion_junctions`, and `exclusion_junctions`. When `object` is a
  `MatisseObject` this defaults to `GetEventData(object)`.

- min_coverage:

  Integer. Minimum total reads per cell per event to report a PSI value.
  Default: `5`.

- na_fill:

  Numeric. Replacement for low-coverage entries. Default: `NA_real_`.

- verbose:

  Logical. Print progress. Default: `TRUE`.

- ...:

  Additional arguments passed to the method.

## Value

A `MatisseObject` (when given one) or a PSI matrix.

- `MatisseObject`: the input object with the `"psi"` assay populated
  inside the embedded Seurat object.

- matrix: a dense matrix (cells × events) of PSI values.

## Details

For each cell \\c\\ and event \\e\\:

\$\$PSI\_{c,e} = \frac{\sum \text{inclusion reads}} {\sum
\text{inclusion reads} + \sum \text{exclusion reads}}\$\$

Results are stored inside the embedded Seurat object as a `Assay5` named
`"psi"`, with:

- `"data"` layer: PSI values in \\\[0,1\]\\ (events × cells).

- `"counts"` layer: inclusion read counts (events × cells).

- `"exclusion"` layer: exclusion read counts (events × cells).

Entries where total coverage falls below `min_coverage` are set to `NA`
in the `"data"` layer.

## See also

[`ComputeIsoformQC`](https://avisrilab.github.io/Matisse/reference/ComputeIsoformQC.md),
[`PlotHeatmap`](https://avisrilab.github.io/Matisse/reference/PlotHeatmap.md)
