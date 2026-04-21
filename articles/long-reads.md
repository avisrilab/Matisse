# Long-read workflow: transcript-based PSI

## When to use this workflow

Use this workflow when your data was generated with a long-read platform
(PacBio, Oxford Nanopore) or when you have used a short-read isoform
quantifier such as **Bagpiper**, **FLAMES**, or **LIQA** that produces
per-transcript per-cell counts. These tools assign each read to a
specific full-length transcript rather than counting at the junction
level.

If you have STARsolo junction counts from standard 10x Chromium short
reads, see the [Short-read
workflow](https://avisrilab.github.io/Matisse/articles/short-reads.md)
instead.

------------------------------------------------------------------------

## What you need

| Input                       | Description                                                                                                                                                                          |
|-----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Seurat object**           | Already processed: gene-expression normalisation, UMAP, cluster labels.                                                                                                              |
| **Transcript count matrix** | *Transcripts × cells* matrix of UMI counts. Row names are transcript IDs matching those in your annotation.                                                                          |
| **SUPPA2 IOE files**        | One or more `.ioe` files from SUPPA2’s `generateEvents` command, one per event type (SE, RI, SS, MX, FL). These map transcripts to inclusion/exclusion sets for each splicing event. |

------------------------------------------------------------------------

## Step 1 — Build the Matisse object

``` r
library(Matisse)

# transcript_counts: transcripts × cells sparse matrix (e.g. from Bagpiper)
# ioe_files: SUPPA2 .ioe output files
obj <- CreateMatisseObjectFromTranscripts(
  seurat            = seu,
  transcript_counts = transcript_counts,
  ioe_files         = c(
    "events_SE.ioe",   # skipped exons
    "events_RI.ioe",   # retained introns
    "events_SS.ioe"    # alternative splice sites
  ),
  min_coverage = 5L
)
```

PSI is calculated automatically during construction — no separate
[`CalculatePSI()`](https://avisrilab.github.io/Matisse/reference/CalculatePSI.md)
call is needed. Each cell × event entry is the fraction of transcripts
carrying the included form of that exon:

$$PSI_{c,e} = \frac{\sum\text{inclusion transcript counts}}{\sum\text{inclusion counts} + \sum\text{exclusion counts}}$$

Cells with fewer than `min_coverage` total transcript counts for a given
event are left as `NA`.

------------------------------------------------------------------------

## Step 2 — Quality control

Check that each cell has enough transcript coverage before moving on to
clustering.

``` r
obj <- ComputeIsoformQC(obj)
PlotQCMetrics(obj, group_by = "cell_type")
```

Remove low-coverage cells and uninformative events:

``` r
obj <- FilterCells(obj, min_pct_covered = 10)
obj <- FilterEvents(obj, min_cells_covered = 20, min_psi_variance = 0.01)
```

------------------------------------------------------------------------

## Step 3 — Normalise transcript counts

For long-read data the transcript-level count matrix often has high
technical variability. `SCTransformTranscripts()` runs SCTransform on
the `"transcript"` assay and then PCA with more components than the
default — useful because isoform variation is subtler than
gene-expression variation and benefits from a wider PCA space.

``` r
# Normalise, regress out sequencing depth, and run PCA in one step.
# Increase n_pca_dims if clusters look under-resolved.
obj <- SCTransformTranscripts(obj, n_pca_dims = 50)
```

------------------------------------------------------------------------

## Step 4 — Cluster and embed

Standard Seurat clustering functions work directly on the Matisse
object. The PSI data and all other slots are preserved throughout.

``` r
obj <- RunUMAP(obj, dims = 1:50)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 0.5)
```

------------------------------------------------------------------------

## Step 5 — Visualise splicing patterns

### Overlay PSI on the UMAP

Each dot is a cell coloured by its PSI value for one event: blue = exon
skipped (low PSI), red = exon included (high PSI).

``` r
PlotPSIUMAP(
  obj,
  event_id = "SE:chr18:3433647-3436055:+",
  title    = "PTBP1 exon 9 — PSI per cell"
)
```

### Compare splicing between cell types

``` r
PlotPSIViolin(
  obj,
  event_id = "SE:chr18:3433647-3436055:+",
  group_by = "cell_type"
)
```

### Survey all events at once

``` r
PlotPSIHeatmap(obj, group_by = "seurat_clusters", max_cells = 400)
```

> **Note on event IDs:** SUPPA2 event IDs use the format
> `TYPE:chr:coords:strand` (e.g. `SE:chr18:100-200:300-400:+`). These
> become the column names in `GetPSI(obj)`, so use the same string when
> calling `PlotPSIUMAP()` or subsetting with `obj[, event_id]`.

------------------------------------------------------------------------

## Accessing data at any point

``` r
# Extract the underlying Seurat object
seu <- GetSeurat(obj)

# Cell metadata via $ (forwarded from Seurat)
obj$seurat_clusters
obj$cell_type

# Full PSI matrix as a sparse matrix (cells × events)
psi <- GetPSI(obj)

# Raw transcript counts (transcripts × cells)
tx <- GetTranscriptCounts(obj)

# Subset to one cell type — PSI and gene expression stay in sync
neurons <- obj[obj$cell_type == "Neuron", ]
```

------------------------------------------------------------------------

## Session info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
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
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.5.3    tools_4.5.3       ragg_1.5.2        bslib_0.10.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.2.0       fs_2.1.0          htmlwidgets_1.6.4
```
