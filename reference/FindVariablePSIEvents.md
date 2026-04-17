# Select highly variable PSI events

Ranks splice events by their cross-cell variability and retains the top
`n_features`. Two ranking methods are available:

- variance:

  Raw per-event variance (after NA imputation). Good default when events
  have comparable mean PSI levels.

- dispersion:

  Coefficient of variation (SD / mean), which down-weights events whose
  variance is high only because their mean is near 0 or 1. Useful when
  events span a wide range of mean PSI.

Selected events are stored in `object@misc\$variable_psi_events` and are
picked up automatically by
[`RunPSIPCA`](https://avisrilab.github.io/Matisse/reference/RunPSIPCA.md)
when `events = NULL`.

## Usage

``` r
FindVariablePSIEvents(
  object,
  n_features = 2000L,
  method = c("variance", "dispersion"),
  verbose = TRUE
)
```

## Arguments

- object:

  A `MatisseObject` with a non-`NULL` `psi` slot.

- n_features:

  Integer. Number of variable events to retain. Default: `2000`.

- method:

  Character. Ranking criterion: `"variance"` (default) or `"dispersion"`
  (coefficient of variation).

- verbose:

  Logical. Default: `TRUE`.

## Value

The input `MatisseObject` with selected event IDs stored in
`object@misc\$variable_psi_events`.

## See also

[`VariablePSIEvents`](https://avisrilab.github.io/Matisse/reference/VariablePSIEvents.md),
[`RunPSIPCA`](https://avisrilab.github.io/Matisse/reference/RunPSIPCA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- FindVariablePSIEvents(obj, n_features = 500)
obj <- RunPSIPCA(obj)   # uses variable events automatically
} # }
```
