library(dplyr)
library(stringr)
library(tidyr)

remove_dataset_indices <- function(dataset){
  
  dataset <- dataset[,-c(1)]
  colnames(dataset) <- dataset[1,]
  dataset <- dataset[-c(1),]
  for(i in 2:ncol(dataset)){
    dataset[,i] <- as.numeric(dataset[,i])
  }
  return(dataset)

}

arrange_samples <- function(dataset, sample_attr, attrs){
  
  sample_attr <- sample_attr[is.element(sample_attr$Label, colnames(dataset)),]
  sample_attr <- dplyr::arrange_at(sample_attr, attrs)
  dataset <- dataset[,c(1, match(sample_attr$Label, colnames(dataset)))]
  return(dataset)
  
}

subset_noise_genes <- function(dataset, skip1, min_percentile){
  
  mean_expr <- rowMeans(dataset[,(skip1+1):ncol(dataset)], na.rm=T)
  subset_vec <- mean_expr>=quantile(mean_expr, min_percentile)
  cat("\nSubset to", sum(subset_vec), "features out of", nrow(dataset), "total features\n")
  return(subset_vec)
  
}

subset_multi_gene_probes <- function(dataset, skip1, max_probes){
  
  top_probe <- dataset %>%
    dplyr::mutate(Mean_Expr=rowMeans(dataset[,(skip1+1):ncol(dataset)], na.rm=T)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::top_n(
      n=max_probes, wt=mean_expr
    )
  
  subset_vec <- is.element(dataset[,c(1)], top_probe$UNIQUE.ID)
  cat("\nSubset to", sum(subset_vec), "features out of", nrow(dataset), "total features\n")
  return(top_probe$idx)
  
}

