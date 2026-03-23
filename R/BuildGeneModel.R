#' Build a gene model from a GTF-derived GRanges object.
#'
#' @param gtf A `GenomicRanges::GRanges` object containing GTF/GFF-style features
#'   with a `type` metadata column and gene/transcript annotations in `mcols(gtf)`.
#' @param gene Character(1). Gene symbol/name to extract from `mcols(gtf)$gene_name`.
#'
#' @return A list containing gene-, transcript-, exon-, and bin-level structures for
#'   the requested gene.
#' @export
BuildGeneModel <- function(
    gtf,
    gene
) {
  
  # ---- Helper checks ----
  .check_granges <- function(x, nm) {
    if (!methods::is(x, "GRanges")) {
      stop(nm, " must be a GenomicRanges::GRanges object.", call. = FALSE)
    }
  }
  
  .check_scalar_character <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1).", call. = FALSE)
    }
  }
  
  # ---- Validate inputs ----
  .check_granges(gtf, "gtf")
  .check_scalar_character(gene, "gene")
  
  gtf.mcols <- S4Vectors::mcols(gtf)
  
  if (!("type" %in% colnames(gtf.mcols))) {
    stop("gtf must contain a 'type' metadata column.", call. = FALSE)
  }
  
  if (!("gene_name" %in% colnames(gtf.mcols))) {
    stop("gtf must contain a 'gene_name' metadata column.", call. = FALSE)
  }
  
  if (!any(gtf.mcols$type == "gene" & gtf.mcols$gene_name == gene, na.rm = TRUE)) {
    stop("No 'gene' feature found in gtf for gene: ", gene, call. = FALSE)
  }
  
  if (!any(gtf.mcols$type == "transcript" & gtf.mcols$gene_name == gene, na.rm = TRUE)) {
    stop("No 'transcript' features found in gtf for gene: ", gene, call. = FALSE)
  }
  
  if (!any(gtf.mcols$type == "exon" & gtf.mcols$gene_name == gene, na.rm = TRUE)) {
    stop("No 'exon' features found in gtf for gene: ", gene, call. = FALSE)
  }
  
  # ---- subset GTF to gene of interest ----
  gene.gr <- gtf[gtf$type == "gene" & mcols(gtf)$gene_name == gene]
  tx.gr <- gtf[gtf$type == "transcript" & mcols(gtf)$gene_name == gene]
  exon.gr <- gtf[gtf$type == "exon" & mcols(gtf)$gene_name == gene]
  
  exon.mcols <- S4Vectors::mcols(exon.gr)
  tx.mcols <- S4Vectors::mcols(tx.gr)
  
  if (!("transcript_id" %in% colnames(exon.mcols))) {
    stop("Exon features for the requested gene must contain a 'transcript_id' metadata column.", call. = FALSE)
  }
  
  if (anyNA(exon.mcols$transcript_id) || any(!nzchar(as.character(exon.mcols$transcript_id)))) {
    stop("Exon features contain missing or empty transcript_id values.", call. = FALSE)
  }
  
  if (!("transcript_id" %in% colnames(tx.mcols))) {
    stop("Transcript features for the requested gene must contain a 'transcript_id' metadata column.", call. = FALSE)
  }
  
  if (anyNA(tx.mcols$transcript_id) || any(!nzchar(as.character(tx.mcols$transcript_id)))) {
    stop("Transcript features contain missing or empty transcript_id values.", call. = FALSE)
  }
  
  # ---- extract transcript identifier fields ----
  tx.id <- mcols(exon.gr)$transcript_id
  tx.name <- mcols(exon.gr)$transcript_name
  if (is.null(tx.name)) {tx.name <- tx.id}
  
  # ---- add transcript identifier fields as columns ----
  mcols(exon.gr)$tx_id <- tx.id
  mcols(exon.gr)$tx_name <- tx.name
  
  # ---- group exons by transcript ----
  exons.by.tx <- split(exon.gr, mcols(exon.gr)$tx_id)
  
  # ---- create non-overlapping genomic bins from all exons ----
  exon.bins <- disjoin(exon.gr, ignore.strand = FALSE)
  
  # ---- annotate each bin with transcript membership ----
  hits <- findOverlaps(exon.bins, exon.gr, ignore.strand = FALSE)
  
  bin.to.tx.id <- split(
    mcols(exon.gr)$tx_id[subjectHits(hits)],
    queryHits(hits)
  )
  
  bin.to.tx.name <- split(
    mcols(exon.gr)$tx_name[subjectHits(hits)],
    queryHits(hits)
  )
  
  # ---- create lists to store bin membership ----
  tx.id.list <- vector("list", length(exon.bins))
  tx.name.list <- vector("list", length(exon.bins))
  
  tx.id.list[as.integer(names(bin.to.tx.id))] <- lapply(bin.to.tx.id, unique)
  tx.name.list[as.integer(names(bin.to.tx.name))] <- lapply(bin.to.tx.name, unique)
  
  tx.id.list <- lapply(tx.id.list, function(x) {
    if (is.null(x)) character(0) else x
  })
  
  tx.name.list <- lapply(tx.name.list, function(x) {
    if (is.null(x)) character(0) else x
  })
  
  # ---- add metadata to bin information ----
  mcols(exon.bins)$gene_name <- gene
  mcols(exon.bins)$bin_id <- paste0(gene, "_bin", seq_along(exon.bins))
  mcols(exon.bins)$tx_ids <- CharacterList(tx.id.list)
  mcols(exon.bins)$tx_names <- CharacterList(tx.name.list)
  mcols(exon.bins)$n_tx <- lengths(mcols(exon.bins)$tx_ids)
  mcols(exon.bins)$is_shared <- mcols(exon.bins)$n_tx > 1
  
  names(exon.bins) <- mcols(exon.bins)$bin_id
  
  # ---- build transcript-to-bin membership mapping table ----
  tx.bin.map <- do.call(
    rbind,
    lapply(seq_along(exon.bins), function(i) {
      txs <- as.character(mcols(exon.bins)$tx_ids[[i]])
      if (length(txs) == 0) return(NULL)
      data.frame(
        tx_id = txs,
        bin_id = mcols(exon.bins)$bin_id[i],
        stringsAsFactors = FALSE
      )
    })
  )
  
  # ---- build transcript metadata table ----
  tx.df <- unique(data.frame(
    tx_id = mcols(tx.gr)$transcript_id,
    tx_name = if ("transcript_name" %in% colnames(mcols(tx.gr))) {
      mcols(tx.gr)$transcript_name
    } else {
      mcols(tx.gr)$transcript_id
    },
    stringsAsFactors = FALSE
  ))
  
  # ---- return ----
  list(
    gene = gene, # gene
    gene_range = gene.gr, # GRanges object of the gene range
    transcripts = tx.gr, # GRanges object of the transcript features
    exons = exon.gr, # GRanges object of the exon features
    exons_by_transcript = exons.by.tx, # GRanges object of each exon split by transcript
    bins = exon.bins, # GRanges object of genomic bins
    transcript_table = tx.df, # Transcript ID to unique transcript name mapping
    transcript_bin_map = tx.bin.map # Transcript ID genomic bin mapping
  )
}