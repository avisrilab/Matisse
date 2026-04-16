# See how your cells splice genes

Matisse measures alternative splicing at single-cell resolution and
layers it on top of your existing Seurat workflow — no re-processing, no
new object format to learn.

[Install](#installation) [View walkthrough
→](https://avisrilab.github.io/Matisse/articles/intro.md)

What you can discover

Biological questions Matisse is built to answer

``` R
<div class="mat-card-num">1</div>
<h3>Cell-type-specific splicing</h3>
<p>Do my neurons and astrocytes process this exon differently?</p>
```

``` R
<div class="mat-card-num">2</div>
<h3>Isoform switches along a trajectory</h3>
<p>Is there a coordinated splicing change as my cells differentiate?</p>
```

``` R
<div class="mat-card-num">3</div>
<h3>Context for bulk RNA-seq</h3>
<p>I see a splicing difference in bulk data — which cell type is driving it?</p>
```

``` R
<div class="mat-card-num">4</div>
<h3>Dataset-wide survey</h3>
<p>Which of the hundreds of alternative exons in my data vary across clusters?</p>
```

![PSI values for Ptbp1 exon 9 overlaid on a UMAP of mouse cortex cells.
Neurons consistently skip this exon; astrocytes include
it.](reference/figures/ptbp1_umap.png)

PSI values for *Ptbp1* exon 9 on a UMAP of mouse cortex cells. Neurons
(left) consistently skip this exon (low PSI); astrocytes (right) include
it (high PSI). Same gene — different splicing. Visible at single-cell
resolution.

Works with your existing setup

Matisse layers on top of Seurat — your clusters, UMAP, and cell labels
stay intact

``` R
<span class="mat-compat-label">Short-read (10x)</span>
<span class="mat-badge">STAR / STARsolo</span>
<span class="mat-badge">junction count matrix</span>
```

``` R
<span class="mat-compat-label">Long-read / isoform</span>
<span class="mat-badge">Bagpiper</span>
<span class="mat-badge">FLAMES</span>
<span class="mat-badge">LIQA</span>
<span class="mat-badge">PacBio MAS-seq</span>
```

``` R
<span class="mat-compat-label">Event annotations</span>
<span class="mat-badge">SUPPA2 generateEvents</span>
<span class="mat-badge">rMATS</span>
<span class="mat-badge">BuildSimpleEvents()</span>
```

Installation

``` r
install.packages("remotes")
remotes::install_github("avisrilab/Matisse")
```

Ready to explore your splicing data?

[View the full walkthrough
→](https://avisrilab.github.io/Matisse/articles/intro.md)
