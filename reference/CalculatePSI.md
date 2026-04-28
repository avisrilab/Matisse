# Calculate PSI matrix from junction or transcript counts

Computes a Percent Spliced In (PSI) matrix for all splice events and
stores it in the `"psi"` assay. Works in both junction mode and
transcript mode. Call this after
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md).

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
  [`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md),
  or a sparse matrix (cells x junctions) when computing PSI outside the
  object.

- events:

  A `data.frame` with columns `event_id`, `inclusion_junctions`, and
  `exclusion_junctions`. Defaults to `GetEventData(object)` (populated
  at construction).

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

- matrix: a dense matrix (cells x events) of PSI values.

## Details

For each cell \\c\\ and event \\e\\:

\$\$PSI\_{c,e} = \frac{\sum \text{inclusion reads}} {\sum
\text{inclusion reads} + \sum \text{exclusion reads}}\$\$

Results are stored inside the embedded Seurat object as `Assay5("psi")`,
with:

- `"data"` layer: PSI values in \\\[0,1\]\\ (events x cells).

- `"counts"` layer: inclusion read counts (events x cells).

- `"exclusion"` layer: exclusion read counts (events x cells).

Entries where total coverage falls below `min_coverage` are set to `NA`
in the `"data"` layer.

`nPercent_isoform` is also written to cell metadata: the percentage of
splice events with a non-`NA` PSI value in each cell.

## See also

[`FilterCells`](https://avisrilab.github.io/Matisse/reference/FilterCells.md),
[`FilterEvents`](https://avisrilab.github.io/Matisse/reference/FilterEvents.md),
[`PlotHeatmap`](https://avisrilab.github.io/Matisse/reference/PlotHeatmap.md)
