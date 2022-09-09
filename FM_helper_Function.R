library(dplyr)
library(stringr)
library(tidyr)

remove_dataset_indices <- function(dataset){
  
  dataset <- dataset[,-c(1)]
  colnames(dataset) <- dataset[1,]
  dataset <- dataset[-c(1),]
  for (i in 2:ncol(dataset)){
    dataset[,i] <- as.numeric(dataset[,i])
  }
  return(dataset)

}

arrange_samples <- function(dataset, sample_attr, attrs){
  
  sample_attr <- sample_attr[sample_attr$Label %in% colnames(dataset),]
  sample_attr <- dplyr::arrange_at(sample_attr, attrs)
  dataset <- dataset[,c(1, match(sample_attr$Label, colnames(dataset)))]
  return(dataset)
  
}

map2SymbolModSnapshot <- function(dataset, dataset_attr, unique_id_col, mapping_tables_dir){

  unique_id <- dataset_attr$Value[dataset_attr$Attribute=="Unique Identifier"]
  if(unique_id=="PROBEID"){
    platform <- dataset_attr$Value[dataset_attr$Attribute=="Mapping Tables"]
    tables_dir <- file.path("Microarray_mapping_tables/rearranged", paste0(platform, "_PROBEID_mapping_tables"))
  } else {
    tables_dir <- file.path("Homo_sapiens_mapping_tables", unique_id)
  }
  mapping_table <- fread(list.files(pattern="SYMBOL", path=file.path(mapping_tables_dir, tables_dir), full.names=T), data.table=F)
  mapping_table <- mapping_table %>% tidyr::separate_rows(SYMBOL, sep=" \\| ")
  mapped_ids <- merge(dataset[,unique_id_col], mapping_table, by=1, all.x=T)
  colnames(mapped_ids)[1] <- "UNIQUE.ID"
  mapped_ids$SYMBOL[is.na(mapped_ids$SYMBOL)] <- mapped_ids$UNIQUE.ID[is.na(mapped_ids$SYMBOL)]
  if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
    mapped_ids <- mapped_ids %>%
      dplyr::arrange(SYMBOL) %>%
      dplyr::group_by(UNIQUE.ID) %>%
      dplyr::slice(1)
  }
  mapped_ids <- mapped_ids[match(dataset[,unique_id_col], mapped_ids$UNIQUE.ID),]
  if(!identical(mapped_ids$UNIQUE.ID, dataset[,unique_id_col])){
    stop("!identical(mapped_ids$UNIQUE.ID, dataset[,unique_id_col])")
  }
  dataset_mapped <- data.frame(mapped_ids, dataset[,2:ncol(dataset)])
  colnames(dataset_mapped)[2] <- "Gene"
  return(dataset_mapped)
  
}

subset_noise_genes <- function(dataset, n_meta_cols, min_percentile){
  
  mean_expr <- rowMeans(dataset[,(n_meta_cols+1):ncol(dataset)], na.rm=T)
  subset_vec <- mean_expr>=quantile(mean_expr, min_percentile)
  cat("\nSubset to", sum(subset_vec), "features out of", nrow(dataset), "total features\n")
  return(subset_vec)
  
}

subset_multi_gene_probes <- function(dataset, n_meta_cols, max_probes){
  
  dataset$mean_expr <- rowMeans(dataset[,(n_meta_cols+1):ncol(dataset)], na.rm=T)
  top_probe <- dataset %>%
    dplyr::mutate(mean_expr=rowMeans(dataset[,(n_meta_cols+1):ncol(dataset)], na.rm=T)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::top_n(n=max_probes, wt=mean_expr)
  
  subset_vec <- dataset[,1] %in% top_probe$UNIQUE.ID
  cat("\nSubset to", sum(subset_vec), "features out of", nrow(dataset), "total features\n")
  return(subset_vec)
  
}

