#' Load single-cell transcript-level RNA-seq data as a sparse matrix.
#'
#' @param barcode.file Character(1). Path to barcode file.
#' @param feature.file Character(1). Path to feature file.
#' @param matrix.file Character(1). Path to Matrix Market file.
#' @param cell.column Integer(1). Column of `barcode.file` to use for cell names (default = 1).
#' @param feats.column Integer(1). Column of `feature.file` to use for feature names (default = 1).
#' @param transpose Logical(1). If TRUE, transpose the matrix after reading (default = FALSE).
#' @param make.unique.features Logical(1). If TRUE, make duplicated feature names unique (default = FALSE).
#'
#' @return A `Matrix::dgCMatrix` with features as rows and cells as columns.
#' @export
LoadTranscriptMatrix <- function(barcode.file,
                     feature.file,
                     matrix.file,
                     cell.column = 1L,
                     feats.column = 1L,
                     transpose = FALSE,
                     make.unique.features = FALSE) {
  
  # ---- Helper functions ----
  .check_path <- function(x, nm) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop(nm, " must be a non-empty character(1) file path.", call. = FALSE)
    }
    if (!file.exists(x)) {
      stop("File not found: ", x, call. = FALSE)
    }
  }
  
  .check_positive_integer <- function(x, nm) {
    if (length(x) != 1L || is.na(x) || x < 1L || x != as.integer(x)) {
      stop(nm, " must be a positive integer(1).", call. = FALSE)
    }
  }
  
  # ---- Verify parameter inputs ----
  .check_path(barcode.file, "barcode.file")
  .check_path(feature.file, "feature.file")
  .check_path(matrix.file, "matrix.file")
  
  .check_positive_integer(cell.column, "cell.column")
  .check_positive_integer(feats.column, "feats.column")
  
  if (!is.logical(transpose) || length(transpose) != 1L || is.na(transpose)) {
    stop("transpose must be a logical(1).", call. = FALSE)
  }
  
  if (!is.logical(make.unique.features) || length(make.unique.features) != 1L || is.na(make.unique.features)) {
    stop("make.unique.features must be a logical(1).", call. = FALSE)
  }
  
  # ---- Read matrix file ----
  mat <- Matrix::readMM(matrix.file)
  
  if (isTRUE(transpose)) {
    mat <- Matrix::t(mat)
  }
  
  # Read barcode & feature file ----
  feats.df <- utils::read.delim(feature.file, header = FALSE, stringsAsFactors = FALSE)
  bcs.df   <- utils::read.delim(barcode.file, header = FALSE, stringsAsFactors = FALSE)
  
  if (nrow(feats.df) < 1L) stop("feature.file has no rows.", call. = FALSE)
  if (nrow(bcs.df) < 1L) stop("barcode.file has no rows.", call. = FALSE)
  
  if (feats.column > ncol(feats.df)) {
    stop("feats.column exceeds number of columns in feature.file.", call. = FALSE)
  }
  if (cell.column > ncol(bcs.df)) {
    stop("cell.column exceeds number of columns in barcode.file.", call. = FALSE)
  }
  
  feats <- trimws(as.character(feats.df[[feats.column]]))
  bcs   <- trimws(as.character(bcs.df[[cell.column]]))
  
  if (anyNA(feats) || any(!nzchar(feats))) {
    stop("feature IDs contain NA or empty values.", call. = FALSE)
  }
  if (anyNA(bcs) || any(!nzchar(bcs))) {
    stop("barcodes contain NA or empty values.", call. = FALSE)
  }
  
  if (nrow(mat) != length(feats)) {
    stop(
      "Row mismatch: matrix has ", nrow(mat), " rows but feature.file has ",
      length(feats), " entries. If cells are rows and features are columns, use transpose = TRUE.",
      call. = FALSE
    )
  }
  
  if (ncol(mat) != length(bcs)) {
    stop(
      "Column mismatch: matrix has ", ncol(mat), " columns but barcode.file has ",
      length(bcs), " entries. If cells are rows and features are columns, use transpose = TRUE.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(bcs)) {
    stop("Duplicate barcodes detected. Cell names must be unique.", call. = FALSE)
  }
  
  if (anyDuplicated(feats)) {
    if (isTRUE(make.unique.features)) {
      warning("Duplicate feature IDs detected; making feature names unique.", call. = FALSE)
      feats <- make.unique(feats)
    } else {
      warning("Duplicate feature IDs detected.", call. = FALSE)
    }
  }
  
  if (!inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }
  
  # ---- Assign RowNames + ColNames and return the matrix ----
  rownames(mat) <- feats
  colnames(mat) <- bcs
  
  mat
}