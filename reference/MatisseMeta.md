# Get or set cell-level isoform metadata

Get or set cell-level isoform metadata

## Usage

``` r
MatisseMeta(object, ...)

MatisseMeta(object) <- value

# S4 method for class 'MatisseObject'
MatisseMeta(object, ...)

# S4 method for class 'MatisseObject'
MatisseMeta(object) <- value
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments (unused).

- value:

  A `data.frame` of cell-level metadata to assign (for the setter).
  Rownames must match cell barcodes.

## Value

For the getter: a `data.frame`. For the setter: the updated
`MatisseObject`.
