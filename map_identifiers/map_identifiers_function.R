library(data.table)
library(tidyr)
library(dplyr)

hs_dir <- "Homo_sapiens_mapping_tables"
array_dir <- "Microarray_mapping_tables/rearranged"

map2Any <- function(
  features, 
  unique_id=c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"), 
  map_to=c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"), 
  unique_id_col, 
  platform=NULL, 
  mapping_tables_dir, 
  keep_all_mappings=F
){
  
  tables_dir <- file.path(hs_dir, unique_id)
  
  if(unique_id=="PROBEID"){
    tables_dir <- file.path(array_dir, paste0(platform, "_PROBEID_mapping_tables"))
  } 
  
  mapping_table <- fread(
    list.files(pattern=map_to, path=file.path(mapping_tables_dir, tables_dir), full.names=T), data.table=F
  )
  mapping_table[,1] <- toupper(mapping_table[,1])
  mapping_table <- data.frame(tidyr::separate_rows(mapping_table, !!as.name(unique_id), sep=" \\| "))
  
  features_orig <- features[,unique_id_col]
  features[,unique_id_col] <- toupper(features[,unique_id_col])
  colnames(features)[colnames(features)=="SYMBOL"] <- "SYMBOL.y" ## Rename any existing SYMBOL column
  
  mapped_ids <- merge(data.frame(features[,unique_id_col]), mapping_table, by.x=1, by.y=1, all.x=T)
  colnames(mapped_ids)[1] <- "UNIQUE.ID"
  
  if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
    
    if(keep_all_mappings==F){
      
      ## If a key maps to multiple identifiers, select only one:
      
      mapped_ids <- mapped_ids %>% 
        dplyr::group_by(UNIQUE.ID) %>% 
        dplyr::slice(1)
      
    } else {
      
      ## Else, collapse rows with one-to-many mappings:
      
      mapped_ids <- mapped_ids %>%
        dplyr::group_by(UNIQUE.ID) %>%
        dplyr::summarise(
          TEMP=paste(!!as.name(map_to), sep=" | ")
        ) %>%
        as.data.frame()
      
      colnames(mapped_ids)[grep("TEMP", colnames(mapped_ids))] <- map_to
      
    }
    
  } ## if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
  
  ## Enforce that data is returned in its original order:
  
  mapped_ids <- mapped_ids[match(features[,unique_id_col], mapped_ids$UNIQUE.ID),]
  
  if(!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)){
    stop("!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)")
  }
  
  features <- data.frame(
    UNIQUE.ID=features_orig, TEMP=mapped_ids[,2], features[,-c(unique_id_col)]
  )
  colnames(features)[2] <- map_to
  
  return(features)
  
}

mapAlias2Symbol <- function(
  features, 
  unique_id_col, 
  mapping_tables_dir, 
  keep_all_mappings=F, 
  fill_NAs=F
){
  
  mapping_table <- fread(
    list.files(pattern="ALIAS", path=file.path(mapping_tables_dir, hs_dir, "SYMBOL"), full.names=T), data.table=F
  )
  mapping_table <- data.frame(separate_rows(mapping_table, ALIAS, sep=" \\| "))
  mapping_table[,2] <- toupper(mapping_table[,2])
  
  features_orig <- features[,unique_id_col]
  features[,unique_id_col] <- toupper(features[,unique_id_col])
  
  mapped_ids <- merge(
    data.frame(features[,unique_id_col]), mapping_table, by.x=1, by.y=2, all.x=T
    )
  
  colnames(mapped_ids)[1] <- "UNIQUE.ID"
  
  if(fill_NAs==T){
    mapped_ids$SYMBOL[is.na(mapped_ids$SYMBOL)] <- mapped_ids$UNIQUE.ID[is.na(mapped_ids$SYMBOL)]
  }
  
  ## For aliases that map to 2+ symbols:
  
  if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
    
    if(keep_all_mappings==F){
      
      ## If a key maps to multiple identifiers, select only one:
      
      mapped_ids <- mapped_ids %>%
        dplyr::group_by(UNIQUE.ID) %>%
        dplyr::slice(1)
      
    } else { ## if(keep_all_mappings==F){
      
      ## Else, collapse rows with one-to-many mappings:
      
      mapped_ids <- mapped_ids %>% 
        dplyr::group_by(UNIQUE.ID) %>%
        dplyr::summarise(SYMBOL=paste(SYMBOL, collapse=" | ")) %>%
        as.data.frame()
      
    } ## if(keep_all_mappings==F){} else {
    
  } ## if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
  
  ## Enforce that data is returned in its original order:
  
  mapped_ids <- mapped_ids[match(features[,unique_id_col], mapped_ids$UNIQUE.ID),]
  
  if(!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)){
    stop("!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)")
  }
  
  features <- data.frame(
    UNIQUE.ID=features_orig, SYMBOL=mapped_ids$SYMBOL, features[,-c(unique_id_col)]
  )
  
  return(features)
  
}


