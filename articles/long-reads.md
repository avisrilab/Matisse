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

| Input | Description |
|----|----|
| **Transcript count matrix** | A 10x-style MatrixMarket triplet (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, optionally `.gz`) from a long-read / isoform quantifier (Bagpiper, FLAMES, LIQA). [`ReadTranscriptMatrix()`](https://avisrilab.github.io/Matisse/reference/ReadTranscriptMatrix.md) loads it in any orientation. |
| **Seurat object** | A Seurat object whose cells match the transcript matrix. For transcript-only data you can build a shell from the matrix and cluster on it directly (Steps 3-4). |
| **SUPPA2 IOE files** | One or more `.ioe` files from SUPPA2’s `generateEvents` command, one per event type (SE, RI, SS, MX, FL). These map transcripts to the splicing events they support. |

------------------------------------------------------------------------

## Step 0 – Read the transcript matrix

[`ReadTranscriptMatrix()`](https://avisrilab.github.io/Matisse/reference/ReadTranscriptMatrix.md)
reads the quantifier’s MatrixMarket triplet (gzip-aware) and returns a
**transcripts × cells** sparse matrix with the dimnames Matisse expects.
Quantifiers disagree on orientation (Bagpiper writes *cells ×
transcripts*; others the reverse) – `orientation = "auto"` detects and
transposes as needed. `strip_version = TRUE` removes a trailing
`.<version>` from transcript IDs (summing any that collapse): SUPPA2
`.ioe` events and the quantifier must use the **same** transcript-ID
namespace, and a silent version mismatch would zero every PSI value.

``` r

library(Matisse)

txc <- ReadTranscriptMatrix(
  "bagpiper/mats",          # dir with matrix.mtx(.gz) + barcodes + features
  cells = colnames(seu)     # optional: subset to your called cells
)
```

------------------------------------------------------------------------

## Step 1 – Build the Matisse object (PSI is computed in the same call)

Pass `transcript_counts` and `events` together. The `events` parameter
accepts either a character vector of SUPPA2 `.ioe` file paths (parsed
internally) or a pre-built `data.frame`. Matisse stores the transcript
counts in an internal assay and computes PSI as part of construction.

``` r

# txc: transcripts x cells matrix from ReadTranscriptMatrix() (Step 0)
# events: SUPPA2 .ioe file paths, one per event type
obj <- CreateMatisseObject(
  seurat            = seu,
  transcript_counts = txc,
  events            = c(
    "events_SE.ioe",   # skipped exons
    "events_RI.ioe",   # retained introns
    "events_SS.ioe"    # alternative splice sites
  ),
  min_coverage      = 5L
)
```

For each cell and event, Matisse sums transcript counts from inclusion
isoforms and exclusion isoforms, then computes:

``` math
PSI_{c,e} = \frac{\sum \text{inclusion transcript counts}}
                   {\sum \text{inclusion counts} + \sum \text{exclusion counts}}
```

Cells with fewer than `min_coverage` total transcripts for a given event
are reported as `NA`. After construction, three QC columns sit in cell
metadata:

- `nCount_isoform` – total transcript counts per cell
- `nFeature_isoform` – number of transcripts with ≥1 read per cell
- `nPercent_isoform` – percentage of splice events with a non-NA PSI per
  cell

To re-compute PSI with different parameters later, call
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
faceted panel. This is the fastest way to spot cells with very few
transcripts or low splicing coverage.

``` r

PlotViolin(obj)
```

### Remove low-quality cells

Remove cells that have too few transcripts or too little splicing
coverage:

``` r

obj <- FilterCells(
  obj,
  min_features_isoform = 5,    # at least 5 transcripts detected
  min_pct_isoform      = 10    # PSI covered for >= 10% of splice events
)
```

### Remove uninformative splice events

Remove events that are observed in very few cells:

``` r

obj <- FilterEvents(obj, min_cells_covered = 20)
```

------------------------------------------------------------------------

## Step 3 – Normalise transcript counts

Normalise the transcript-level count matrix before clustering.
SCTransform applies variance stabilisation and corrects for sequencing
depth.

``` r

obj <- SCTransform(obj)
obj <- RunPCA(obj, assay = "SCT", npcs = 50)
```

------------------------------------------------------------------------

## Step 4 – Cluster and embed

Standard Seurat clustering functions work directly on the Matisse
object. The PSI data and all other slots are preserved throughout.

``` r

obj <- RunUMAP(obj, dims = 1:50)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 0.5)
```

------------------------------------------------------------------------

## Step 5 – Visualise splicing patterns

### Overlay PSI on the UMAP

Each dot is a cell coloured by its PSI value for one splice event: blue
= exon skipped (low PSI), red = exon included (high PSI).

``` r

PlotUMAP(
  obj,
  feature = "SE:chr18:3433647-3436055:+",
  title   = "PTBP1 exon 9 -- PSI per cell"
)
```

### Compare splicing between cell types

``` r

PlotViolin(
  obj,
  feature  = "SE:chr18:3433647-3436055:+",
  group_by = "cell_type"
)
```

### Survey all events at once

``` r

PlotHeatmap(obj, group_by = "seurat_clusters", max_cells = 400)
```

### Inspect isoform usage with a sashimi plot

[`PlotSashimi()`](https://avisrilab.github.io/Matisse/reference/PlotSashimi.md)
draws junction arcs scaled by aggregate read count, coloured by role
(inclusion = blue, exclusion = red). Facet by cell type to compare
isoform usage across populations.

``` r

PlotSashimi(
  obj,
  event_id = "SE:chr18:3433648-3434699:3434801-3436055:-",
  group_by = "cell_type",
  title    = "PTBP1 exon 9 -- coverage by cell type"
)
```

> **Note on event IDs:** SUPPA2 SE event IDs use the format
> `SE:chr:donor1-acceptor1:donor2-acceptor2:strand`
> (e.g. `SE:chr18:3433648-3434699:3434801-3436055:-`). These become the
> feature names in `GetPSI(obj)`, so use the same string when calling
> [`PlotUMAP()`](https://avisrilab.github.io/Matisse/reference/PlotUMAP.md),
> [`PlotSashimi()`](https://avisrilab.github.io/Matisse/reference/PlotSashimi.md),
> or subsetting.

------------------------------------------------------------------------

## Accessing data at any point

``` r

# Extract the underlying Seurat object
seu <- GetSeurat(obj)

# Cell metadata via $ (forwarded from Seurat)
obj$seurat_clusters
obj$cell_type

# Full PSI matrix as a sparse matrix (cells x splice events)
psi <- GetPSI(obj)

# Raw transcript counts (transcripts x cells) via the native Seurat layer API
tx <- SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]], "counts")

# Subset to one cell type -- PSI and gene expression stay in sync
neurons <- obj[obj$cell_type == "Neuron", ]
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
