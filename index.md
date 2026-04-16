# Understand your cells, layer by layer

Matisse brings isoform-resolved splicing and chromatin accessibility
into your single-cell workflow — on top of Seurat and Signac, using the
same cells, the same clusters, the same UMAP.

[Install](#installation) [View walkthrough
→](https://avisrilab.github.io/Matisse/articles/intro.md)

What you can discover

Questions Matisse is built to answer

1

### Cell-type-specific splicing

Do my neurons and astrocytes process this exon differently — and by how
much?

2

### Chromatin shapes isoforms

Is the splicing switch I see linked to chromatin accessibility changes
at the same locus?

3

### Isoform switches along a trajectory

Is there a coordinated splicing change as my cells differentiate or
respond to a stimulus?

4

### Context for bulk RNA-seq

I see a splicing difference in bulk data — which cell type is
responsible?

![PSI values for Ptbp1 exon 9 overlaid on a UMAP of mouse cortex cells.
Neurons skip this exon; astrocytes include
it.](reference/figures/ptbp1_umap.png)

PSI values for *Ptbp1* exon 9 on a UMAP of mouse cortex cells. Neurons
(blue) consistently skip this exon (low PSI); astrocytes (red) include
it (high PSI). Same gene — different isoforms — visible at single-cell
resolution.

Works with your existing setup

Matisse layers on top of Seurat and Signac — your clusters, UMAP, and
cell labels stay intact

Short-read RNA (10x) STAR / STARsolo junction count matrix

Long-read / isoform Bagpiper FLAMES LIQA PacBio MAS-seq

ATAC / chromatin 10x Multiome Signac ArchR

Event annotations SUPPA2 generateEvents rMATS BuildSimpleEvents()

Installation

``` r
install.packages("remotes")
remotes::install_github("avisrilab/Matisse")
```

Ready to explore your data?

[View the full walkthrough
→](https://avisrilab.github.io/Matisse/articles/intro.md)
