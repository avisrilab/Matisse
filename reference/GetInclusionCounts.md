# Get inclusion read count matrix

Retrieves inclusion counts from the `"counts"` layer of the `"psi"`
`ChromatinAssay`.

## Usage

``` r
GetInclusionCounts(object, ...)

# S4 method for class 'MatisseObject'
GetInclusionCounts(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments (unused).

## Value

A sparse matrix (cells × events) of inclusion read counts.
