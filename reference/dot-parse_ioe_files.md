# Parse one or more SUPPA2 IOE files into an event table

Parse one or more SUPPA2 IOE files into an event table

## Usage

``` r
.parse_ioe_files(ioe_files)
```

## Arguments

- ioe_files:

  Character vector of file paths.

## Value

A data.frame with columns: `event_id`, `gene_id`, `chr`, `strand`,
`event_type`, `inclusion_transcripts`, `exclusion_transcripts`.
