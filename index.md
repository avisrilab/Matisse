---
pagetitle: "Matisse — See how your cells splice genes"
---

<!-- badges: start -->
[![R-CMD-check](https://github.com/avisrilab/Matisse/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/avisrilab/Matisse/actions/workflows/R-CMD-check.yml)
[![pkgdown](https://github.com/avisrilab/Matisse/actions/workflows/pkgdown.yml/badge.svg)](https://avisrilab.github.io/Matisse)
<!-- badges: end -->

<div class="mat-hero">
<h1>See how your cells<br>splice genes</h1>
<p class="mat-hero-sub">
  Matisse measures alternative splicing at single-cell resolution and
  layers it on top of your existing Seurat workflow —
  no re-processing, no new object format to learn.
</p>
<a href="#installation" class="btn mat-btn-primary">Install</a>
<a href="articles/intro.html" class="btn mat-btn-outline">View walkthrough &rarr;</a>
</div>

<section class="mat-section">
<p class="mat-section-title">What you can discover</p>
<p class="mat-section-sub">Biological questions Matisse is built to answer</p>

<div class="mat-cards">
  <div class="mat-card">
    <div class="mat-card-num">1</div>
    <h3>Cell-type-specific splicing</h3>
    <p>Do my neurons and astrocytes process this exon differently?</p>
  </div>
  <div class="mat-card">
    <div class="mat-card-num">2</div>
    <h3>Isoform switches along a trajectory</h3>
    <p>Is there a coordinated splicing change as my cells differentiate?</p>
  </div>
  <div class="mat-card">
    <div class="mat-card-num">3</div>
    <h3>Context for bulk RNA-seq</h3>
    <p>I see a splicing difference in bulk data — which cell type is driving it?</p>
  </div>
  <div class="mat-card">
    <div class="mat-card-num">4</div>
    <h3>Dataset-wide survey</h3>
    <p>Which of the hundreds of alternative exons in my data vary across clusters?</p>
  </div>
</div>
</section>

<figure class="mat-figure">
<img src="man/figures/ptbp1_umap.png"
     alt="PSI values for Ptbp1 exon 9 overlaid on a UMAP of mouse cortex cells. Neurons consistently skip this exon; astrocytes include it.">
<figcaption>
  PSI values for <em>Ptbp1</em> exon 9 on a UMAP of mouse cortex cells.
  Neurons (left) consistently skip this exon (low PSI); astrocytes (right) include it (high PSI).
  Same gene — different splicing. Visible at single-cell resolution.
</figcaption>
</figure>

<section class="mat-section">
<p class="mat-section-title">Works with your existing setup</p>
<p class="mat-section-sub">Matisse layers on top of Seurat — your clusters, UMAP, and cell labels stay intact</p>

<div class="mat-compat">
  <div class="mat-compat-group">
    <span class="mat-compat-label">Short-read (10x)</span>
    <span class="mat-badge">STAR / STARsolo</span>
    <span class="mat-badge">junction count matrix</span>
  </div>
  <div class="mat-compat-group">
    <span class="mat-compat-label">Long-read / isoform</span>
    <span class="mat-badge">Bagpiper</span>
    <span class="mat-badge">FLAMES</span>
    <span class="mat-badge">LIQA</span>
    <span class="mat-badge">PacBio MAS-seq</span>
  </div>
  <div class="mat-compat-group">
    <span class="mat-compat-label">Event annotations</span>
    <span class="mat-badge">SUPPA2 generateEvents</span>
    <span class="mat-badge">rMATS</span>
    <span class="mat-badge">BuildSimpleEvents()</span>
  </div>
</div>
</section>

<section class="mat-section" id="installation">
<p class="mat-section-title">Installation</p>

<div class="mat-install">

```r
install.packages("remotes")
remotes::install_github("avisrilab/Matisse")
```

</div>
</section>

<div class="mat-cta-footer">
  <p>Ready to explore your splicing data?</p>
  <a href="articles/intro.html" class="btn mat-btn-primary">View the full walkthrough &rarr;</a>
</div>
