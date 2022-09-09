## Ref: https://stackoverflow.com/questions/24614391/intersect-all-possible-combinations-of-list-elements

source("/home/rebecca/code/misc/jaccard_index.R")

overlap_matrix <- function(list, jaccard=F){
  
  if(jaccard){
    
    jac_idx <- lapply(list, function(x) unlist(
      lapply(list, function(y){
        jaccard_index(x, y)*100
      })
    ))
    
    mat <- do.call(cbind, jac_idx)
    
  } else {
    
    mat <- crossprod(table(stack(list)))
    
  }
  
  return(as.matrix(mat))
}


# ## Divide all elements by the union of all elements:
# union_genes <- n_distinct(do.call(c, dataset_genes))
# mat <- mat/union_genes*100