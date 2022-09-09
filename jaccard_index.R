jaccard_index <- function(x, y) {
  
  intrsct <- length(intersect(x, y))
  union <- length(x) + length(y) - intrsct
  return(intrsct/union)
  
}