# Build a junction-ID event table from SUPPA2 events

SUPPA2 `.ioe` files describe each splicing event with *transcript* IDs,
but junction mode needs the inclusion/exclusion features to be
*junction* IDs that match a STARsolo SJ matrix. The genomic coordinates
of every defining junction are, however, fully encoded in the SUPPA2
`event_id` (e.g. `SE:chr1:804222-804776:804966-807217:+`). This function
parses that grammar per event type and emits the `chr-start-end-strand`
junction IDs Matisse expects, returning the `data.frame` that
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
accepts via `events`. This is a deterministic coordinate adapter, not a
heuristic event caller.

## Usage

``` r
BuildJunctionEvents(
  ioe_files,
  event_types = c("SE", "A3", "A5", "MX", "AF", "AL"),
  junction_universe = NULL,
  intron_offset = "auto",
  verbose = TRUE
)
```

## Arguments

- ioe_files:

  Character vector of SUPPA2 `.ioe` file path(s).

- event_types:

  Character vector of event types to keep. Default:
  `c("SE", "A3", "A5", "MX", "AF", "AL")`.

- junction_universe:

  Optional character vector of observed junction IDs (e.g.
  `colnames(ReadSTARsoloSJ(...))`). Used to auto-calibrate
  `intron_offset` and to report match rate. Default: `NULL`.

- intron_offset:

  Either `"auto"` (calibrate against `junction_universe`; falls back to
  `c(1L, 1L)` if absent) or an integer vector `c(left, right)` added to
  the donor / subtracted from the acceptor exon-boundary coordinate.
  Default: `"auto"`.

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

A `data.frame` with columns `event_id`, `gene_id`, `chr`, `strand`,
`event_type`, `inclusion_features`, `exclusion_features` (the last two
semicolon-separated junction IDs).

## Details

**Coordinate convention.** SUPPA2 reports exon-boundary positions (last
base of the donor exon, first base of the acceptor exon). STARsolo SJ
reports the first/last *intronic* base. The two differ by a fixed offset
(intron = `[donor + 1, acceptor - 1]`). When `junction_universe` is
supplied and `intron_offset = "auto"`, the offset is calibrated by
maximising overlap with the observed STARsolo junctions; otherwise the
standard `c(1L, 1L)` is used.

**Supported event types.** `SE`, `A3`, `A5`, `MX`, `AF`, `AL`. `RI`
(retained intron) is *not* junction-quantifiable - intron retention has
no supporting split junction - and is dropped with a warning if
requested.

## See also

[`ReadSTARsoloSJ`](https://avisrilab.github.io/Matisse/reference/ReadSTARsoloSJ.md),
[`CreateMatisseObject`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md).

## Examples

``` r
if (FALSE) { # \dontrun{
jxn    <- ReadSTARsoloSJ("Solo.out/SJ/raw", cells = colnames(seu))
events <- BuildJunctionEvents(
  list.files("events", "\\.ioe$", full.names = TRUE),
  junction_universe = colnames(jxn)
)
obj <- CreateMatisseObject(seu, junction_counts = jxn, events = events)
} # }
```
