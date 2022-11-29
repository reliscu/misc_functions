## Ref: https://stackoverflow.com/questions/24614391/intersect-all-possible-combinations-of-list-elements

source("/home/rebecca/code/misc/jaccard_index.R")

overlap_matrix <- function(genes_list, jaccard=F){
  
  if(jaccard){
    
    jac_idx <- lapply(genes_list, function(x) unlist(lapply(genes_list, function(y){jaccard_index(x, y)*100})))
    mat <- do.call(cbind, jac_idx)
    
  } else {
    mat <- crossprod(table(stack(genes_list)))
  }
  
  return(as.matrix(mat))
  
}
