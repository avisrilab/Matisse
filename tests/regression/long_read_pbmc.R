# ===========================================================================
# Regression: long-read (Bagpiper transcript) pipeline on a PIP-seq PBMC
# dataset. Exercises ReadTranscriptMatrix -> CreateMatisseObject (transcript
# mode) -> PSI -> QC -> all five plots, with hard assertions.
#
#   Rscript reg_testing/long_reads/long_read_pbmc.R          # uses cache
#   FRESH=1 Rscript reg_testing/long_reads/long_read_pbmc.R  # rebuild cache
#
# Clones the short-read scaffold (reg_testing/short_reads/short_read_pbmc.R)
# and the long-reads.Rmd clustering sequence. Not part of the package build
# (reg_testing/ is gitignored).
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
pkgload::load_all(".", quiet = TRUE)
set.seed(1)

root     <- "reg_testing/long_reads"
mats_dir <- file.path(root, "mats")
ev_dir   <- "reg_testing/short_reads/events"      # reuse SUPPA2 .ioe
fig_dir  <- file.path(root, "figs")
cache    <- file.path(root, ".cache_inputs.rds")
dir.create(fig_dir, showWarnings = FALSE)

say <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), sprintf(...)))

fresh <- nzchar(Sys.getenv("FRESH")) || !file.exists(cache)

if (!fresh) {
  say("Loading cached called-cell transcript matrix (%s) ...", cache)
  txc <- readRDS(cache)$txc
} else {
  # --- 1. Read Bagpiper transcript matrix (transcripts x cells) ------------
  say("Reading Bagpiper transcript matrix (this is large) ...")
  txc_full <- ReadTranscriptMatrix(mats_dir)        # auto-orients + strips ver
  say("Full matrix: %d transcripts x %d barcodes",
      nrow(txc_full), ncol(txc_full))

  # --- 2. Cell calling via UMI knee (Kneedle on the rank plot) -------------
  tot <- Matrix::colSums(txc_full)
  ord <- order(tot, decreasing = TRUE)
  pos <- tot[ord][tot[ord] > 0]
  x   <- log10(seq_along(pos)); y <- log10(as.numeric(pos))
  xn  <- (x - min(x)) / (max(x) - min(x))
  yn  <- (y - min(y)) / (max(y) - min(y))
  knee <- which.max(abs(yn - (1 - xn)))             # max dist to anti-diagonal
  knee <- max(500L, min(knee, 60000L))              # sane PBMC-scale clamp
  called <- names(pos)[seq_len(knee)]
  say("UMI knee at rank %d (min UMIs=%d); calling %d cells",
      knee, as.integer(pos[knee]), length(called))

  txc <- txc_full[, called, drop = FALSE]
  txc <- txc[Matrix::rowSums(txc) > 0, , drop = FALSE]   # drop empty tx
  saveRDS(list(txc = txc), cache)
  say("Cached called-cell matrix (%d tx x %d cells) -> %s",
      nrow(txc), ncol(txc), cache)
}

# --- 3. Seurat shell + Matisse object (transcript mode) --------------------
ioe <- list.files(ev_dir, pattern = "\\.ioe$", full.names = TRUE)
stopifnot(length(ioe) > 0L)

# Transcript-ID overlap diagnostic (analog of short-read calibration):
# fraction of .ioe transcripts observed in the Bagpiper feature space.
ioe_tx <- unique(unlist(lapply(ioe, function(f) {
  d <- utils::read.delim(f, sep = "\t", stringsAsFactors = FALSE,
                         quote = "", comment.char = "")
  unlist(strsplit(d$total_transcripts, ",", fixed = TRUE))
}), use.names = FALSE))
ioe_tx     <- sub("\\.[0-9]+$", "", ioe_tx)
tx_overlap <- mean(ioe_tx %in% rownames(txc))
say("IOE-transcript overlap with Bagpiper features: %.1f%% (%d/%d)",
    100 * tx_overlap, sum(ioe_tx %in% rownames(txc)), length(ioe_tx))

seu <- CreateSeuratObject(txc, project = "pbmc_lr", min.cells = 1)
say("Seurat shell: %d transcripts x %d cells", nrow(seu), ncol(seu))

obj <- CreateMatisseObject(seu, transcript_counts = txc,
                           events = ioe, min_coverage = 5L)

# --- 4. Cluster + embed: reuse the long-reads.Rmd dispatch sequence --------
obj <- SCTransform(obj, verbose = FALSE)            # transcript mode -> SCT
obj <- RunPCA(obj, assay = "SCT", npcs = 50, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:50, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
say("Clusters: %d", length(unique(GetSeurat(obj)$seurat_clusters)))

# --- 5. Coarse PBMC cell typing at transcript resolution -------------------
# No tx2gene matrix (user choice: cluster on transcripts directly); build
# ENST->ENSG from the .ioe gene_id/total_transcripts columns, then score
# canonical marker ENSGs by their summed-transcript expression per cluster.
tx2gene <- {
  m <- do.call(rbind, lapply(ioe, function(f) {
    d <- utils::read.delim(f, sep = "\t", stringsAsFactors = FALSE,
                           quote = "", comment.char = "")
    g <- sub(";.*", "", d$gene_id)
    tx <- strsplit(d$total_transcripts, ",", fixed = TRUE)
    data.frame(tx = sub("\\.[0-9]+$", "", unlist(tx)),
               gene = rep(g, lengths(tx)), stringsAsFactors = FALSE)
  }))
  setNames(m$gene[!duplicated(m$tx)], m$tx[!duplicated(m$tx)])
}
marker_ensg <- list(
  `T`      = "ENSG00000167286",  # CD3D
  B        = "ENSG00000156738",  # MS4A1
  Mono     = c("ENSG00000170458", "ENSG00000090382"),  # CD14, LYZ
  `NK`     = c("ENSG00000105374", "ENSG00000115523"),  # NKG7, GNLY
  DC       = "ENSG00000179639",  # FCER1A
  Platelet = "ENSG00000163736"   # PPBP
)
tx_counts <- SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]],
                                        layer = "counts")  # tx x cells
libsize   <- Matrix::colSums(tx_counts)
clusters  <- GetSeurat(obj)$seurat_clusters
score <- sapply(marker_ensg, function(gs) {
  tx <- names(tx2gene)[tx2gene %in% gs]
  tx <- intersect(tx, rownames(tx_counts))
  if (length(tx) == 0L) return(rep(0, length(unique(clusters))))
  cpm <- log1p(Matrix::colSums(tx_counts[tx, , drop = FALSE]) /
               pmax(libsize, 1) * 1e4)
  tapply(cpm, clusters, mean)[levels(factor(clusters))]
})
if (is.list(score)) score <- do.call(cbind, score)
rownames(score) <- levels(factor(clusters))
score <- scale(score)                                  # z across clusters
score[is.na(score)] <- 0
cl2type <- colnames(score)[max.col(score, ties.method = "first")]
names(cl2type) <- rownames(score)
# MatisseObject has no `$<-`; the supported path is AddMetaData (delegates
# to Seurat, aligning by barcode rownames).
ct_df <- data.frame(
  cell_type = unname(cl2type[as.character(clusters)]),
  row.names = colnames(GetSeurat(obj)))
obj <- AddMetaData(obj, ct_df)
say("Cell types: %s",
    paste(sprintf("%s=%d",
                  names(table(GetSeurat(obj)$cell_type)),
                  as.integer(table(GetSeurat(obj)$cell_type))),
          collapse = ", "))

# --- 6. RAW PSI stats (slot-level; reuse short-read idiom) -----------------
psi      <- GetPSI(obj)                          # cells x events sparse
xs       <- psi@x
not_na   <- !is.na(xs)
n_non_na <- sum(not_na)
vals     <- xs[not_na]
frac_var <- if (length(vals)) mean(vals > 0 & vals < 1) else 0
col_of_x      <- rep.int(seq_len(ncol(psi)), diff(psi@p))
cov_per_event <- as.integer(
  tapply(not_na, factor(col_of_x, levels = seq_len(ncol(psi))), sum))
cov_per_event[is.na(cov_per_event)] <- 0L
names(cov_per_event) <- colnames(psi)
n_informative <- sum(cov_per_event >= 20)
say("RAW PSI: %d cells x %d events | non-NA=%s | in (0,1)=%.1f%% | events w/ >=20 covered=%d",
    nrow(psi), ncol(psi), format(n_non_na, big.mark = ","),
    100 * frac_var, n_informative)

psi_col <- function(m, ev) {
  j <- match(ev, colnames(m)); out <- rep(NA_real_, nrow(m))
  rng <- (m@p[j] + 1L):m@p[j + 1L]
  if (m@p[j + 1L] > m@p[j]) out[m@i[rng] + 1L] <- m@x[rng]
  out
}

# --- 7. Biological signal: cell-type-differential splicing -----------------
ct       <- GetSeurat(obj)$cell_type
inf_ev   <- utils::head(
  names(sort(cov_per_event[cov_per_event >= 20], decreasing = TRUE)), 200L)
ct_range <- vapply(inf_ev, function(e) {
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

# --- 8. Focus on informative events, then plot -----------------------------
obj <- FilterEvents(obj, min_cells_covered = 20)
obj <- FilterCells(obj, min_features_isoform = 5)
psi_f <- GetPSI(obj)
say("After filtering: %d cells x %d events", nrow(psi_f), ncol(psi_f))

cf_col <- rep.int(seq_len(ncol(psi_f)), diff(psi_f@p))
cov_f  <- as.integer(tapply(!is.na(psi_f@x),
            factor(cf_col, levels = seq_len(ncol(psi_f))), sum))
cov_f[is.na(cov_f)] <- 0L
names(cov_f) <- colnames(psi_f)
top_event <- names(sort(cov_f, decreasing = TRUE))[1]
say("Top-covered event: %s (%d covered cells)", top_event, max(cov_f))
# PlotSashimi (transcript mode) supports only SE / RI event types; pick the
# top-covered SE/RI event for that panel specifically.
sash_pool     <- cov_f[grepl("^(SE|RI):", names(cov_f))]
sashimi_event <- if (length(sash_pool))
  names(sort(sash_pool, decreasing = TRUE))[1] else top_event
say("Sashimi event (SE/RI): %s", sashimi_event)

ggsave2 <- function(f, p, ...)
  ggplot2::ggsave(file.path(fig_dir, f), p, ...)
ok_plot <- function(tag, expr) {
  p <- tryCatch(expr, error = function(e) {
    say("PLOT FAIL [%s]: %s", tag, conditionMessage(e)); NULL })
  if (!is.null(p)) ggsave2(paste0(tag, ".png"), p, width = 8, height = 5)
  !is.null(p)
}
plots_ok <- c(
  qc      = ok_plot("qc_violin",   PlotViolin(obj)),
  umap    = ok_plot("umap_event",  PlotUMAP(obj, feature = top_event)),
  violin  = ok_plot("violin_event",
                     PlotViolin(obj, feature = top_event,
                                group_by = "cell_type")),
  heatmap = ok_plot("heatmap",
                     PlotHeatmap(obj, group_by = "cell_type",
                                 max_cells = 400)),
  sashimi = ok_plot("sashimi",
                     PlotSashimi(obj, event_id = sashimi_event,
                                 group_by = "cell_type"))
)
say("Figures (%d/5 rendered) -> %s", sum(plots_ok), fig_dir)

# --- 9. Hard assertions (success criteria) ---------------------------------
errs <- character(0)
chk  <- function(cond, msg) if (!isTRUE(cond)) errs <<- c(errs, msg)

chk(tx_overlap >= 0.30,
    sprintf("IOE/Bagpiper transcript overlap %.1f%% < 30%% (namespace mismatch)",
            100 * tx_overlap))
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
saveRDS(obj, file.path(root, "pbmc_lr_matisse.rds"))
