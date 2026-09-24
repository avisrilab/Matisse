---
pagetitle: "Matisse: per-cell splicing (PSI) in Seurat"
---

<div class="mat-hero">
<h1>Understand your cells,<br>layer by layer</h1>
<p class="mat-hero-sub">Matisse computes percent spliced in (PSI) for each splicing event in each cell and stores it as an assay in your Seurat object, next to gene expression: the same cells, the same clusters, the same UMAP.</p>
<a href="#installation" class="btn mat-btn-primary">Install</a>&nbsp;<a href="articles/intro.html" class="btn mat-btn-outline">View walkthrough &rarr;</a>
</div>

<div class="mat-section">
<p class="mat-section-title">What you can discover</p>
<p class="mat-section-sub">Questions Matisse is built to answer</p>
<div class="mat-cards">
<div class="mat-card">
<div class="mat-card-num">1</div>
<h3>Cell-type-specific splicing</h3>
<p>Do two cell types include this exon at different rates, and by how much?</p>
</div>
<div class="mat-card">
<div class="mat-card-num">2</div>
<h3>Long reads or short reads</h3>
<p>Start from per-cell transcript counts from a long-read quantifier such as Bagpiper, or from STARsolo splice-junction counts from 10x short reads.</p>
</div>
</div>
</div>

<div class="mat-section">
<p class="mat-section-title">Works with your existing setup</p>
<p class="mat-section-sub">Matisse layers on top of Seurat: your clusters, UMAP, and cell labels stay intact</p>
<div class="mat-compat">
<div class="mat-compat-group">
<span class="mat-compat-label">Short-read RNA (10x)</span>
<span class="mat-badge">STAR / STARsolo</span>
<span class="mat-badge">junction count matrix</span>
</div>
<div class="mat-compat-group">
<span class="mat-compat-label">Long-read / isoform</span>
<span class="mat-badge">Bagpiper</span>
<span class="mat-badge">FLAMES</span>
<span class="mat-badge">LIQA</span>
</div>
<div class="mat-compat-group">
<span class="mat-compat-label">Event annotations</span>
<span class="mat-badge">SUPPA2 generateEvents</span>
</div>
</div>
</div>

<div class="mat-section" id="installation">
<p class="mat-section-title">Installation</p>
<div class="mat-install">
<pre><code class="language-r">install.packages("remotes")
remotes::install_github("avisrilab/Matisse")</code></pre>
</div>
</div>

<div class="mat-cta-footer">
<p>Ready to explore your data?</p>
<a href="articles/intro.html" class="btn mat-btn-primary">View the full walkthrough &rarr;</a>
</div>
