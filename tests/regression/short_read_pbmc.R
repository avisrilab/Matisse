# ===========================================================================
# Regression: short-read (STARsolo junction) pipeline on a 10x PBMC dataset.
# Exercises ReadSTARsoloSJ -> BuildJunctionEvents -> CreateMatisseObject ->
# PSI -> QC -> all four plots, with hard assertions on correctness.
#
#   Rscript reg_testing/short_reads/short_read_pbmc.R          # uses cache
#   FRESH=1 Rscript reg_testing/short_reads/short_read_pbmc.R  # rebuild cache
#
# Not part of the package build (reg_testing/ is gitignored).
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
pkgload::load_all(".", quiet = TRUE)
set.seed(1)

root     <- "reg_testing/short_reads"
solo     <- file.path(root, "star_Solo.out")
gene_dir <- file.path(solo, "Gene", "filtered")
sj_dir   <- file.path(solo, "SJ", "raw")
ev_dir   <- file.path(root, "events")
fig_dir  <- file.path(root, "figs")
cache    <- file.path(root, ".cache_inputs.rds")
dir.create(fig_dir, showWarnings = FALSE)

say <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), sprintf(...)))

fresh <- nzchar(Sys.getenv("FRESH")) || !file.exists(cache)

if (!fresh) {
  say("Loading cached Seurat + SJ inputs (%s) ...", cache)
  ci  <- readRDS(cache)
  seu <- ci$seu
  jxn <- ci$jxn
} else {
  # --- 1. Base Seurat object from the gene matrix --------------------------
  say("Reading Gene/filtered ...")
  gene <- ReadMtx(
    mtx      = file.path(gene_dir, "matrix.mtx"),
    cells    = file.path(gene_dir, "barcodes.tsv"),
    features = file.path(gene_dir, "features.tsv"),
    feature.column = 2L
  )
  seu <- CreateSeuratObject(gene, project = "pbmc_sr", min.cells = 3,
                            min.features = 200)
  say("Seurat: %d cells x %d genes", ncol(seu), nrow(seu))

  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

  # --- 2. Coarse PBMC cell typing (canonical markers, per-cluster argmax) --
  marker_sets <- list(
    `T`      = c("CD3D", "CD3E", "IL7R"),
    B        = c("MS4A1", "CD79A", "CD79B"),
    Mono     = c("CD14", "LYZ", "FCN1"),
    `NK`     = c("NKG7", "GNLY", "KLRD1"),
    DC       = c("FCER1A", "CST3"),
    Platelet = c("PPBP", "PF4")
  )
  avg <- AverageExpression(seu, group.by = "seurat_clusters",
                           assays = "RNA", layer = "data",
                           verbose = FALSE)$RNA
  score <- sapply(marker_sets, function(g) {
    g <- intersect(g, rownames(avg))
    if (length(g) == 0) return(rep(0, ncol(avg)))
    rowMeans(scale(t(avg[g, , drop = FALSE])), na.rm = TRUE)
  })
  if (is.list(score)) score <- do.call(cbind, score)
  rownames(score) <- colnames(avg)
  cl2type <- colnames(score)[max.col(score, ties.method = "first")]
  names(cl2type) <- sub("^g", "", rownames(score))
  seu$cell_type <- factor(unname(cl2type[as.character(seu$seurat_clusters)]))
  say("Cell types: %s",
      paste(sprintf("%s=%d", names(table(seu$cell_type)),
                    as.integer(table(seu$cell_type))), collapse = ", "))

  # --- 3. STARsolo junctions ----------------------------------------------
  jxn <- ReadSTARsoloSJ(sj_dir, cells = colnames(seu))
  stopifnot(nrow(jxn) > 0L, ncol(jxn) > 1000L)

  saveRDS(list(seu = seu, jxn = jxn), cache)
  say("Cached inputs -> %s", cache)
}

# --- 4. SUPPA2 events -> junction-ID event table ---------------------------
ioe <- list.files(ev_dir, pattern = "\\.ioe$", full.names = TRUE)
stopifnot(length(ioe) > 0L)
events <- BuildJunctionEvents(ioe, junction_universe = colnames(jxn))

# --- 5. Matisse object + PSI (report BEFORE filtering) ---------------------
obj <- CreateMatisseObject(seu, junction_counts = jxn,
                           events = events, min_coverage = 5L)

# GetPSI -> sparse dgCMatrix, cells x events, NAs stored explicitly in @x.
# Operate on slots only; is.na() on the full matrix would densify (~20 GB).
psi      <- GetPSI(obj)
xs       <- psi@x
not_na   <- !is.na(xs)
n_non_na <- sum(not_na)
vals     <- xs[not_na]
frac_var <- if (length(vals)) mean(vals > 0 & vals < 1) else 0
# covered (= non-NA) entries per event: events are columns -> use @p spans
col_of_x      <- rep.int(seq_len(ncol(psi)), diff(psi@p))
cov_per_event <- as.integer(
  tapply(not_na, factor(col_of_x, levels = seq_len(ncol(psi))), sum))
cov_per_event[is.na(cov_per_event)] <- 0L
names(cov_per_event) <- colnames(psi)
n_informative <- sum(cov_per_event >= 20)
say("RAW PSI: %d cells x %d events | non-NA entries=%s | in (0,1)=%.1f%% | events w/ >=20 covered cells=%d",
    nrow(psi), ncol(psi), format(n_non_na, big.mark = ","),
    100 * frac_var, n_informative)

# Extract one event column as a length-ncell vector (NA where not stored).
psi_col <- function(m, ev) {
  j <- match(ev, colnames(m)); out <- rep(NA_real_, nrow(m))
  rng <- (m@p[j] + 1L):m@p[j + 1L]
  if (m@p[j + 1L] > m@p[j]) out[m@i[rng] + 1L] <- m@x[rng]
  out
}

# --- 6. Biological signal: cell-type-differential splicing -----------------
# Short-read per-cell SJ PSI is intrinsically sparse and bimodal (most
# covered events resolve to 0 or 1 at this depth), so we do NOT gate on a
# fixed fraction of intermediate PSI. Instead we verify the pipeline yields
# real, cell-type-structured splicing: among well-covered events, the mean
# PSI per cell type must vary for at least one event.
ct        <- GetSeurat(obj)$cell_type
inf_ev    <- names(sort(cov_per_event[cov_per_event >= 20],
                        decreasing = TRUE))
inf_ev    <- utils::head(inf_ev, 200L)
ct_range  <- vapply(inf_ev, function(e) {
  mu <- tapply(psi_col(psi, e), ct, function(z) mean(z, na.rm = TRUE))
  mu <- mu[is.finite(mu)]
  if (length(mu) < 2L) 0 else diff(range(mu))
}, numeric(1))
max_ct_range <- if (length(ct_range)) max(ct_range) else 0
best_var_ev  <- if (length(ct_range)) inf_ev[which.max(ct_range)] else NA
if (!is.na(best_var_ev)) {
  mu <- tapply(psi_col(psi, best_var_ev), ct,
               function(z) mean(z, na.rm = TRUE))
  say("Most cell-type-variable event: %s (range=%.2f) | %s",
      best_var_ev, max_ct_range,
      paste(sprintf("%s=%.2f", names(mu), as.numeric(mu)),
            collapse = ", "))
}

# Informational PTPRC (CD45) probe (not a gate; often low-coverage at 10x).
ptprc_ev <- {
  edf <- GetSeurat(obj)[["psi"]][[]]
  e   <- rownames(edf)[edf$gene_id == "ENSG00000081237"]
  e[e %in% colnames(psi)]
}
if (length(ptprc_ev) && max(cov_per_event[ptprc_ev]) > 0) {
  e  <- ptprc_ev[which.max(cov_per_event[ptprc_ev])]
  mu <- tapply(psi_col(psi, e), ct, function(z) mean(z, na.rm = TRUE))
  say("PTPRC %s (cov=%d) | %s", e, cov_per_event[e],
      paste(sprintf("%s=%.2f", names(mu), as.numeric(mu)),
            collapse = ", "))
} else {
  say("PTPRC: no covered event (expected at shallow 10x SJ depth).")
}

# --- 7. Focus on informative events, then plot -----------------------------
obj <- FilterEvents(obj, min_cells_covered = 20)
obj <- FilterCells(obj, min_features_isoform = 5)  # NOT min_pct: 283k
                                                    # genome-wide events make
                                                    # a % threshold meaningless
psi_f <- GetPSI(obj)                          # cells x events (filtered)
say("After filtering: %d cells x %d events", nrow(psi_f), ncol(psi_f))

cf_col  <- rep.int(seq_len(ncol(psi_f)), diff(psi_f@p))
cov_f   <- as.integer(tapply(!is.na(psi_f@x),
             factor(cf_col, levels = seq_len(ncol(psi_f))), sum))
cov_f[is.na(cov_f)] <- 0L
names(cov_f) <- colnames(psi_f)
top_event <- names(sort(cov_f, decreasing = TRUE))[1]
say("Top-covered event: %s (%d covered cells)",
    top_event, max(cov_f))

ggsave2 <- function(f, p, ...)
  ggplot2::ggsave(file.path(fig_dir, f), p, ...)
ok_plot <- function(tag, expr) {
  p <- tryCatch(expr, error = function(e) {
    say("PLOT FAIL [%s]: %s", tag, conditionMessage(e)); NULL })
  if (!is.null(p)) ggsave2(paste0(tag, ".png"), p, width = 8, height = 5)
  !is.null(p)
}
plots_ok <- c(
  qc       = ok_plot("qc_violin",   PlotViolin(obj)),
  umap     = ok_plot("umap_event",  PlotUMAP(obj, feature = top_event)),
  violin   = ok_plot("violin_event",
                      PlotViolin(obj, feature = top_event,
                                 group_by = "cell_type")),
  heatmap  = ok_plot("heatmap",
                      PlotHeatmap(obj, group_by = "cell_type",
                                  max_cells = 400)),
  sashimi  = ok_plot("sashimi",
                      PlotSashimi(obj, event_id = top_event,
                                  group_by = "cell_type"))
)
say("Figures (%d/5 rendered) -> %s", sum(plots_ok), fig_dir)

# --- 8. Hard assertions (success criteria) ---------------------------------
errs <- character(0)
chk  <- function(cond, msg) if (!isTRUE(cond)) errs <<- c(errs, msg)

chk(n_non_na > 50000,
    sprintf("RAW PSI nearly empty (%d non-NA) -> ID-mismatch regression",
            n_non_na))
chk(n_informative >= 100,
    sprintf("only %d events have >=20 covered cells", n_informative))
chk(max_ct_range >= 0.10,
    sprintf("no informative event varies by cell type (max range %.2f)",
            max_ct_range))
chk(all(plots_ok),
    sprintf("plots failed: %s",
            paste(names(plots_ok)[!plots_ok], collapse = ", ")))

if (length(errs) > 0) {
  say("FAIL:\n - %s", paste(errs, collapse = "\n - "))
  quit(status = 1)
}
say("ALL CHECKS PASSED")
saveRDS(obj, file.path(root, "pbmc_sr_matisse.rds"))
