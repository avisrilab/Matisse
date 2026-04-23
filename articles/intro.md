# Getting started with Matisse

## What does Matisse do?

Genes can be spliced in different ways – certain exons may be included
in some transcripts but skipped in others. This is called **alternative
splicing**, and it means one gene can produce multiple protein variants
with different functions.

Matisse measures alternative splicing **one cell at a time**. For each
cell it calculates a **PSI value** (Percent Spliced In) for each
splicing event:

- **PSI = 1** – every transcript in that cell includes the exon
- **PSI = 0** – every transcript skips the exon
- **PSI = 0.5** – half include, half skip

By comparing PSI values across cell types, clusters, or conditions you
can find which populations splice genes differently – and by how much.
Matisse sits on top of your existing
[Seurat](https://satijalab.org/seurat/) workflow. Your gene expression
data, UMAP, and cluster labels stay intact; splicing information is
added alongside them.

------------------------------------------------------------------------

## Two entry points

Matisse supports two data types via a single constructor. What you pass
in determines the mode:

| Data type                                                          | How to create                                               | When to use                                                                 |
|--------------------------------------------------------------------|-------------------------------------------------------------|-----------------------------------------------------------------------------|
| **Long-read / transcript quantification** (Bagpiper, FLAMES, LIQA) | `CreateMatisseObject(transcript_counts=..., ioe_files=...)` | You have a transcripts x cells count matrix and SUPPA2 IOE files            |
| **Short-read 10x Chromium** (STARsolo junction counts)             | `CreateMatisseObject(junction_counts=..., event_data=...)`  | You have a cells x junctions count matrix from STARsolo `--soloFeatures SJ` |

Once the object is built, all downstream functions – QC, filtering,
normalisation, visualisation – are identical for both data types.

------------------------------------------------------------------------

## Quick-start (long reads)

``` r
library(Matisse)

obj <- CreateMatisseObject(
  seurat            = seu,
  transcript_counts = tx_counts,    # transcripts x cells from Bagpiper/FLAMES
  ioe_files         = "events_SE.ioe",
  min_coverage      = 5L
)
# PSI is computed at construction; summarise before clustering
psi_summary <- SummarizePSI(obj)
obj <- ComputeIsoformQC(obj)
# Normalise transcript counts, then reduce and cluster
obj <- SCTransform(obj)
obj <- RunPCA(obj, assay = "SCT", npcs = 50)
obj <- RunUMAP(obj, dims = 1:50)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 0.5)
PlotUMAP(obj, feature = "SE:chr18:3433647-3436055:+")
```

-\> **Full walkthrough:** [Long-read
workflow](https://avisrilab.github.io/Matisse/articles/long-reads.md)

------------------------------------------------------------------------

## Quick-start (short reads)

``` r
library(Matisse)

obj <- CreateMatisseObject(
  seurat          = seu,          # your existing Seurat object
  junction_counts = jxn_counts,   # cells x junctions from STARsolo
  event_data      = event_df      # splice event annotation
)
obj <- CalculatePSI(obj, min_coverage = 5)
obj <- ComputeIsoformQC(obj)
PlotUMAP(obj, feature = "PTBP1:SE:chr18:3433647-3436055")
```

-\> **Full walkthrough:** [Short-read
workflow](https://avisrilab.github.io/Matisse/articles/short-reads.md)

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
