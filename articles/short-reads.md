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

| Input                      | Description                                                                                                                                                                                                                                                            |
|----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Seurat object**          | Already processed: normalisation, UMAP, cluster labels.                                                                                                                                                                                                                |
| **Junction count matrix**  | *Cells × junctions* sparse matrix of read counts. Row names are cell barcodes; column names are junction IDs in the format `chr:start:end`. Produced by STARsolo with `--soloFeatures SJ`.                                                                             |
| **Event annotation table** | A data frame with one row per splicing event, specifying which junctions support inclusion versus exclusion. Can be built with [`BuildSimpleEvents()`](https://avisrilab.github.io/Matisse/reference/BuildSimpleEvents.md) or loaded from a SUPPA2 / rMATS annotation. |

------------------------------------------------------------------------

## Step 1 — Build the Matisse object

``` r
library(Matisse)

obj <- CreateMatisseObject(
  seurat          = seu,
  junction_counts = jxn_counts,   # cells × junctions matrix from STARsolo
  event_data      = event_df,     # splice event annotation table
  junction_data   = junction_df   # optional: genomic coordinates per junction
)
```

If you don’t have a hand-curated event annotation,
[`BuildSimpleEvents()`](https://avisrilab.github.io/Matisse/reference/BuildSimpleEvents.md)
can generate one automatically from the junction table — each junction
becomes its own “simple” event with all other junctions for the same
gene treated as exclusion evidence:

``` r
event_df <- BuildSimpleEvents(junction_df)
```

------------------------------------------------------------------------

## Step 2 — Calculate PSI

For each cell and each event, Matisse sums the reads from
inclusion-supporting junctions and exclusion-supporting junctions, then
computes:

$$PSI_{c,e} = \frac{\sum\text{inclusion junction reads}}{\sum\text{inclusion reads} + \sum\text{exclusion reads}}$$

Cells with fewer than `min_coverage` total reads for an event are left
as `NA`.

``` r
obj <- CalculatePSI(obj, min_coverage = 5)
```

------------------------------------------------------------------------

## Step 3 — Quality control

``` r
obj <- ComputeIsoformQC(obj)
PlotQCMetrics(obj, group_by = "cell_type")
```

Remove cells that are too sparse or events that are covered in too few
cells:

``` r
obj <- FilterCells(
  obj,
  min_junctions      = 5,
  min_junction_reads = 20,
  min_pct_covered    = 10
)

obj <- FilterEvents(
  obj,
  min_cells_covered = 20,
  min_psi_variance  = 0.01
)
```

------------------------------------------------------------------------

## Step 4 — Visualise splicing patterns

### Where does the splicing switch happen on the UMAP?

``` r
PlotPSIUMAP(
  obj,
  event_id = "PTBP1:SE:chr18:3433647-3436055",
  title    = "PTBP1 exon 9 — PSI per cell"
)
```

### Distributions by cell type

``` r
PlotPSIViolin(
  obj,
  event_id = "PTBP1:SE:chr18:3433647-3436055",
  group_by = "cell_type"
)
```

### Overview of all events

``` r
PlotPSIHeatmap(obj, group_by = "cell_type", max_cells = 400)
```

### Inspect raw junction coverage for a gene

``` r
PlotJunctionCoverage(obj, gene_id = "PTBP1", group_by = "cell_type")
```

------------------------------------------------------------------------

## Accessing data at any point

``` r
# Extract the Seurat object for any standard Seurat step
seu <- GetSeurat(obj)

# Cell metadata (forwarded from the embedded Seurat object)
obj$seurat_clusters
obj$cell_type

# PSI matrix: cells × events
psi <- GetPSI(obj)

# Raw junction counts: cells × junctions
jxn <- GetJunctionCounts(obj)

# Subset to a specific cell type
neurons <- obj[obj$cell_type == "Neuron", ]
```

------------------------------------------------------------------------

## Combining datasets

If you have multiple samples and want to analyse them together, merge
the individual Matisse objects before running QC:

``` r
# Process each sample independently up to CreateMatisseObject + CalculatePSI,
# then merge
combined <- MergeMatisse(
  obj_sample1, obj_sample2,
  add_cell_ids = c("S1", "S2")
)
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
