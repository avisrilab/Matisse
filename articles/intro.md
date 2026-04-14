# Introduction to Matisse

## Overview

**Matisse** (Multi-modal Analysis of Transcript Isoforms in Single-cell
Sequencing Experiments) provides an integrated framework for analysing
isoform-resolved single-cell RNA-seq data. It builds on
[Seurat](https://satijalab.org/seurat/) and
[Signac](https://stuartlab.org/signac/), adding:

- A `MatisseObject` that wraps a `Seurat` object and co-stores isoform
  layers
- PSI (Percent Spliced In) matrix calculation from junction read counts
- Isoform-aware quality control metrics and filtering
- Visualization functions for PSI values on UMAPs, violins, and heatmaps

## Installation

``` r
# Install development version from GitHub
# install.packages("remotes")
remotes::install_github("avisrilab/Matisse")
```

## The MatisseObject

The central data structure is a `MatisseObject`. It consists of:

| Slot               | Contents                                                            |
|--------------------|---------------------------------------------------------------------|
| `seurat`           | An embedded `Seurat` object (gene expression, reductions, metadata) |
| `junction_counts`  | Sparse matrix: cells × junctions                                    |
| `psi`              | Sparse matrix: cells × splice events (PSI in \[0, 1\])              |
| `inclusion_counts` | Sparse matrix: cells × events                                       |
| `exclusion_counts` | Sparse matrix: cells × events                                       |
| `event_data`       | Splice event annotation table                                       |
| `junction_data`    | Junction annotation table                                           |
| `isoform_metadata` | Per-cell isoform QC metrics                                         |

## Quick Start

### 1. Prepare your data

Matisse expects:

1.  A `Seurat` object (processed gene-level data from any pipeline)
2.  A sparse junction count matrix (cells × junctions), e.g. from
    STARsolo, Cell Ranger, or a long-read quantifier
3.  A splice event annotation table (e.g. exported from rMATS or SUPPA2)

``` r
library(Matisse)
library(Seurat)
library(Matrix)

# -- Toy data for illustration -----------------------------------------------
set.seed(42)
n_cells <- 200; n_genes <- 500; n_jxns <- 50

gene_counts <- matrix(rpois(n_genes * n_cells, 10),
                      nrow     = n_genes,
                      dimnames = list(paste0("Gene", seq_len(n_genes)),
                                      paste0("CELL", seq_len(n_cells))))

seu <- CreateSeuratObject(gene_counts) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA() |>
  RunUMAP(dims = 1:10) |>
  FindNeighbors(dims = 1:10) |>
  FindClusters()

# Junction counts (cells × junctions)
jxn_counts <- Matrix(
  matrix(rpois(n_cells * n_jxns, lambda = 3),
         nrow     = n_cells,
         dimnames = list(colnames(seu),
                         paste0("chr1:", seq(1000, by = 500, length.out = n_jxns),
                                "-", seq(1500, by = 500, length.out = n_jxns), ":+"))),
  sparse = TRUE
)
```

### 2. Define splice events

``` r
# Build a simple one-versus-rest event table from the junction list
jd <- data.frame(
  junction_id = colnames(jxn_counts),
  chr         = "chr1",
  start       = seq(1000, by = 500, length.out = n_jxns),
  end         = seq(1500, by = 500, length.out = n_jxns),
  strand      = "+",
  gene_id     = rep(paste0("Gene", 1:10), each = 5)
)

# Use the built-in helper for simple (one-vs-rest) events
ed <- BuildSimpleEvents(jd)
head(ed)
```

### 3. Create a MatisseObject

``` r
obj <- CreateMatisseObject(
  seurat          = seu,
  junction_counts = jxn_counts,
  event_data      = ed,
  junction_data   = jd
)
obj
```

### 4. Calculate PSI

``` r
obj <- CalculatePSI(obj, min_coverage = 5)
```

### 5. Quality control

``` r
obj <- ComputeIsoformQC(obj)
PlotQCMetrics(obj, group_by = "seurat_clusters")
```

Filter out cells with insufficient isoform coverage:

``` r
obj <- FilterCells(
  obj,
  min_junctions      = 5,
  min_junction_reads = 20,
  min_pct_covered    = 10
)
```

Remove poorly covered events:

``` r
obj <- FilterEvents(obj, min_cells_covered = 20, min_psi_variance = 0.01)
```

### 6. Visualisation

``` r
# UMAP coloured by PSI of a specific event
PlotPSIUMAP(obj, event_id = "Gene1:chr1:1000-1500:+:chr1:1000-1500:+")

# Violin plot across clusters
PlotPSIViolin(obj,
              event_id = "Gene1:chr1:1000-1500:+:chr1:1000-1500:+",
              group_by = "seurat_clusters")

# Heatmap of top variable events
PlotPSIHeatmap(obj, max_cells = 100)

# Per-gene junction coverage
PlotJunctionCoverage(obj, gene = "Gene1")
```

## Accessing the embedded Seurat object

All Seurat functionality remains directly accessible:

``` r
# Extract the Seurat object
seu_embedded <- GetSeurat(obj)

# Seurat-style metadata access via $
obj$seurat_clusters

# Subset by cell barcode (keeps both Seurat and isoform layers in sync)
subset_obj <- obj[c("CELL1", "CELL2", "CELL3"), ]
```

## Next steps

- **Vignette: PSI analysis** — differential splicing between cell types
- **Vignette: Integration with Signac** — linking splicing to chromatin
  accessibility in multiome data

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
#> [25] rlang_1.2.0       fs_2.0.1          htmlwidgets_1.6.4
```
