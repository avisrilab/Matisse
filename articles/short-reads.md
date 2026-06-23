# Short-read workflow: junction-based PSI

## When to use this workflow

Use this workflow when your data comes from standard 10x Chromium short
reads processed with **STARsolo** or **STAR**. These tools can output a
per-cell junction count table alongside the gene expression matrix.
Matisse uses those junction counts to compute PSI (Percent Spliced In)
for each splicing event in each cell.

If you have full-transcript counts from a long-read platform or a
transcript quantifier (Bagpiper, FLAMES, LIQA), see the [Long-read
workflow](https://avisrilab.github.io/Matisse/articles/long-reads.md)
instead.

------------------------------------------------------------------------

## What you need

| Input | Description |
|----|----|
| **Seurat object** | Already processed: normalisation, UMAP, cluster labels. |
| **STARsolo `SJ` output** | The `Solo.out/SJ/raw` directory from STARsolo run with `--soloFeatures Gene SJ`. Only a *raw* SJ matrix is produced, so it is subset to your filtered cells on read. |
| **SUPPA2 `.ioe` events** | One file per event type from `suppa.py generateEvents` on your annotation GTF. The junction coordinates are encoded in each `event_id`; Matisse converts them to junction IDs for you. Matisse does not ship a heuristic event caller — bring your own annotation. |

------------------------------------------------------------------------

## Step 1 – Load STARsolo junctions and build the event table

[`ReadSTARsoloSJ()`](https://avisrilab.github.io/Matisse/reference/ReadSTARsoloSJ.md)
reads the SJ matrix and relabels each junction into the
`chr-start-end-strand` ID form Matisse uses (so sashimi plots work for
free). Pass your filtered cell barcodes to subset the raw matrix.

``` r

library(Matisse)

jxn <- ReadSTARsoloSJ(
  "Solo.out/SJ/raw",
  cells = colnames(seu)        # filtered barcodes (e.g. from Gene/filtered)
)
```

SUPPA2 `.ioe` files describe events with *transcript* IDs, but junction
mode needs *junction* IDs.
[`BuildJunctionEvents()`](https://avisrilab.github.io/Matisse/reference/BuildJunctionEvents.md)
parses the genomic coordinates encoded in each SUPPA2 `event_id` and
emits the matching junction IDs. It is a deterministic coordinate
adapter, not an event caller. Supported event types: `SE`, `A3`, `A5`,
`MX`, `AF`, `AL`. (`RI` is dropped — intron retention has no supporting
split junction, so it cannot be scored from SJ counts.)

``` r

events <- BuildJunctionEvents(
  list.files("events", pattern = "\\.ioe$", full.names = TRUE),
  junction_universe = colnames(jxn)   # auto-calibrates the exon/intron
)                                     # coordinate offset to your data
```

SUPPA2 reports exon-boundary positions while STARsolo reports intronic
coordinates; passing `junction_universe` lets
[`BuildJunctionEvents()`](https://avisrilab.github.io/Matisse/reference/BuildJunctionEvents.md)
calibrate the fixed offset by maximising overlap with the junctions
actually observed, and reports the match rate so you can sanity-check
the annotation matches your genome build.

## Step 2 – Build the Matisse object (PSI is computed in the same call)

``` r

obj <- CreateMatisseObject(
  seurat          = seu,
  junction_counts = jxn,          # cells x junctions from ReadSTARsoloSJ()
  events          = events,       # junction-ID table from BuildJunctionEvents()
  min_coverage    = 5L            # minimum reads per cell per event for PSI
)
```

[`CreateMatisseObject()`](https://avisrilab.github.io/Matisse/reference/CreateMatisseObject.md)
calculates PSI as part of construction. After it returns, the object
holds three things ready to use:

- `nCount_isoform` – total junction read counts per cell (in cell
  metadata)
- `nFeature_isoform` – number of junctions with ≥1 read per cell (in
  cell metadata)
- `nPercent_isoform` – percentage of splice events with a non-NA PSI
  value per cell (in cell metadata)
- A `Assay5("psi")` holding the PSI matrix and the per-event annotation
  as feature metadata

For each cell and event, Matisse sums the inclusion-supporting junctions
and exclusion-supporting junctions, then computes:

``` math
PSI_{c,e} = \frac{\sum \text{inclusion junction reads}}
                   {\sum \text{inclusion reads} + \sum \text{exclusion reads}}
```

Cells with fewer than `min_coverage` total reads for an event are left
as `NA`. To re-compute PSI with different parameters later, call
`CalculatePSI(obj, min_coverage = ...)`. To skip the PSI step at
construction (rare), pass `defer_psi = TRUE` and call
[`CalculatePSI()`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
manually.

------------------------------------------------------------------------

## Step 2 – Quality control

### Visualise QC metrics

Call
[`PlotViolin()`](https://avisrilab.github.io/Matisse/reference/PlotViolin.md)
with no `feature` argument to automatically plot all three isoform QC
metrics (`nCount_isoform`, `nFeature_isoform`, `nPercent_isoform`) as a
faceted panel.

``` r

PlotViolin(obj)
```

### Remove low-quality cells

``` r

obj <- FilterCells(
  obj,
  min_features_isoform = 5,    # at least 5 junctions detected
  min_pct_isoform      = 10    # PSI covered for >= 10% of splice events
)
```

### Remove uninformative splice events

``` r

obj <- FilterEvents(obj, min_cells_covered = 20)
```

------------------------------------------------------------------------

## Step 3 – Visualise splicing patterns

### Where does the splicing switch happen on the UMAP?

``` r

PlotUMAP(
  obj,
  feature = "SE:chr18:3433648-3434699:3434801-3436055:-",
  title   = "PTBP1 exon 9 -- PSI per cell"
)
```

### Distributions by cell type

``` r

PlotViolin(
  obj,
  feature  = "SE:chr18:3433648-3434699:3434801-3436055:-",
  group_by = "cell_type"
)
```

### Overview of all events

``` r

PlotHeatmap(obj, group_by = "cell_type", max_cells = 400)
```

### Sashimi plot for a specific event

[`PlotSashimi()`](https://avisrilab.github.io/Matisse/reference/PlotSashimi.md)
shows junction arcs scaled by read count, coloured by inclusion (blue)
vs exclusion (red). Use `group_by` to compare cell types side by side.

``` r

PlotSashimi(
  obj,
  event_id  = "SE:chr18:3433648-3434699:3434801-3436055:-",
  group_by  = "cell_type",
  arc_scale = "sqrt",
  title     = "PTBP1 exon 9 -- junction coverage by cell type"
)
```

------------------------------------------------------------------------

## Accessing data at any point

``` r

# Extract the Seurat object for any standard Seurat step
seu <- GetSeurat(obj)

# Cell metadata (forwarded from the embedded Seurat object)
obj$seurat_clusters
obj$cell_type

# PSI matrix: cells x splice events
psi <- GetPSI(obj)

# Raw junction counts: cells x junctions (via the native Seurat layer API)
jxn <- Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]],
                                            "counts"))

# Subset to a specific cell type
neurons <- obj[obj$cell_type == "Neuron", ]
```

------------------------------------------------------------------------

## Combining datasets

If you have multiple samples and want to analyse them together, merge
the individual Matisse objects before running QC:

``` r

# Build each sample's MatisseObject (PSI is computed at construction), then
# merge. Both objects must share the same `input.mode` and event set.
combined <- MergeMatisse(
  obj_sample1, obj_sample2,
  add_cell_ids = c("S1", "S2")
)
```

------------------------------------------------------------------------

## Session info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.59         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.0    tools_4.6.0       ragg_1.5.2        bslib_0.11.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.2.0       fs_2.1.0          htmlwidgets_1.6.4
```
