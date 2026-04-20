# Get transcript count matrix

Retrieves raw transcript counts from the `"transcript"` `Assay5` stored
inside the embedded Seurat object.

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

A sparse matrix (transcripts × cells) of raw counts, or `NULL` if no
`"transcript"` assay exists.
