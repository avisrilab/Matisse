# Add columns to the isoform metadata

Add columns to the isoform metadata

## Usage

``` r
AddIsoformMetadata(object, metadata, ...)

# S4 method for class 'MatisseObject'
AddIsoformMetadata(object, metadata, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- metadata:

  A named `data.frame` or named numeric/character vector. Rownames (or
  names) must match cell barcodes.

- ...:

  Additional arguments (unused).

## Value

The updated `MatisseObject`.

## Methods (by class)

- `AddIsoformMetadata(MatisseObject)`: Add or update columns in the
  isoform metadata.
