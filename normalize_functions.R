library(Matrix)

normalize_fxn <- function(expr, scale_factor=1e6) {
  expr <- sweep(expr,2,colSums(expr),"/")*scale_factor
  return(expr)
}

normalize_sparse <- function(A, scale_factor=1e6) {
  A@x <- A@x/Matrix::colSums(A)[A@j+1L]*scale_factor
  return(A)
  # https://stackoverflow.com/questions/39284774/column-rescaling-for-a-very-large-sparse-matrix-in-r
}

log2_sparse <- function(A) {
  A@x <- log2(A@x+1)
  return(A)
}