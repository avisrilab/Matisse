# Get exclusion read count matrix

Retrieves exclusion counts from the `"exclusion"` layer of the `"psi"`
`ChromatinAssay`.

## Usage

``` r
GetExclusionCounts(object, ...)

# S4 method for class 'MatisseObject'
GetExclusionCounts(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments (unused).

## Value

A sparse matrix (cells × events) of exclusion read counts.
