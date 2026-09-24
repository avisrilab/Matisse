# Understand your cells, layer by layer

Matisse computes percent spliced in (PSI) for each splicing event in
each cell and stores it as an assay in your Seurat object, next to gene
expression: the same cells, the same clusters, the same UMAP.

[Install](#installation) [View walkthrough
→](https://avisrilab.org/Matisse/articles/intro.md)

What you can discover

Questions Matisse is built to answer

1

### Cell-type-specific splicing

Do two cell types include this exon at different rates, and by how much?

2

### Long reads or short reads

Start from per-cell transcript counts from a long-read quantifier such
as Bagpiper, or from STARsolo splice-junction counts from 10x short
reads.

Works with your existing setup

Matisse layers on top of Seurat: your clusters, UMAP, and cell labels
stay intact

Short-read RNA (10x) STAR / STARsolo junction count matrix

Long-read / isoform Bagpiper FLAMES LIQA

Event annotations SUPPA2 generateEvents

Installation

``` r
install.packages("remotes")
remotes::install_github("avisrilab/Matisse")
```

Ready to explore your data?

[View the full walkthrough
→](https://avisrilab.org/Matisse/articles/intro.md)
