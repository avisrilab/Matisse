# Get transcript count matrix (transcript mode only)

Retrieves raw transcript counts from the `"isoform"` `Assay5` in
transcript mode.

## Usage

``` r
GetTranscriptCounts(object, ...)

# S4 method for class 'MatisseObject'
GetTranscriptCounts(object, ...)
```

## Arguments

- object:

  A `MatisseObject`.

- ...:

  Additional arguments (unused).

## Value

A sparse matrix (transcripts x cells) of raw counts, or `NULL` if the
object is in junction mode.
