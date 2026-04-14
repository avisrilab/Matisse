# Calculate PSI matrix from junction counts

Computes a Percent Spliced In (PSI) matrix for all splice events defined
in `event_data`. For each cell \\c\\ and event \\e\\:

## Usage

``` r
CalculatePSI(object, events = NULL, ...)

# S4 method for class 'MatisseObject'
CalculatePSI(
  object,
  events = NULL,
  min_coverage = 5L,
  na_fill = NA_real_,
  verbose = TRUE
)

# S4 method for class 'ANY'
CalculatePSI(
  object,
  events,
  min_coverage = 5L,
  na_fill = NA_real_,
  verbose = TRUE
)
```

## Arguments

- object:

  A
  [`MatisseObject`](https://avisrilab.github.io/Matisse/reference/MatisseObject-class.md)
  that has a non-`NULL` `junction_counts` slot, or a sparse matrix
  (cells x junctions).

- events:

  When `object` is a matrix: a `data.frame` with columns `event_id`,
  `inclusion_junctions`, and `exclusion_junctions`. When `object` is a
  `MatisseObject` this defaults to `GetEventData(object)`.

- ...:

  Additional arguments passed to the method.

- min_coverage:

  Integer. Minimum total reads per cell per event to report a PSI value.
  Default: `5`.

- na_fill:

  Numeric. Replacement for low-coverage entries. Default: `NA_real_`.

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

A `MatisseObject` (when given one) or a PSI matrix.

- `MatisseObject`: the input object with `psi`, `inclusion_counts`, and
  `exclusion_counts` slots populated.

- matrix: a dense matrix (cells x events) of PSI values.

## Details

\$\$PSI\_{c,e} = \frac{\sum \text{inclusion reads}} {\sum
\text{inclusion reads} + \sum \text{exclusion reads}}\$\$

Entries where the total coverage (inclusion + exclusion) falls below
`min_coverage` are set to `na_fill` (default `NA`).

## See also

[`ComputeIsoformQC`](https://avisrilab.github.io/Matisse/reference/ComputeIsoformQC.md),
[`PlotPSIHeatmap`](https://avisrilab.github.io/Matisse/reference/PlotPSIHeatmap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Build a small toy junction matrix
jxn_mat <- Matrix::sparseMatrix(
  i = c(1,1,2,2),
  j = c(1,2,2,3),
  x = c(10, 5, 8, 3),
  dims = c(3, 4),
  dimnames = list(
    paste0("Cell", 1:3),
    c("jxn1","jxn2","jxn3","jxn4")
  )
)
events <- data.frame(
  event_id             = "SE_gene1",
  gene_id              = "gene1",
  chr                  = "chr1",
  strand               = "+",
  event_type           = "SE",
  inclusion_junctions  = "jxn1;jxn2",
  exclusion_junctions  = "jxn3",
  stringsAsFactors     = FALSE
)
psi_mat <- CalculatePSI(jxn_mat, events, min_coverage = 3)
} # }
```
