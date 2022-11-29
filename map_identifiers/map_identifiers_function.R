library(tidyr)
library(dplyr)
library(data.table)

map2Any <- function(features, unique_id=c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"), map_to, unique_id_col, platform=NULL, tables_path="/home/rebecca/omicon/mapping_tables/2022-11-02", keep_all=F){
  
  if(unique_id=="PROBEID"){
    tables_path <- file.path(tables_path, "microarray_mapping_tables", paste0(platform, "_PROBEID_mapping_tables"))
  } else {
    tables_path <- file.path(tables_path, "homo_sapiens_mapping_tables", unique_id)
  }
  
  mapping_table <- fread(list.files(pattern=map_to, path=tables_path, full.names=T), data.table=F)
  mapping_table[,c(1)] <- toupper(mapping_table[,c(1)])
  mapping_table <- data.frame(
    tidyr::separate_rows(mapping_table, !!as.name(unique_id), sep=" \\| ")
  )
  
  features1 <- features[,unique_id_col]
  features[,unique_id_col] <- toupper(features[,unique_id_col])
  colnames(features)[colnames(features)=="SYMBOL"] <- "SYMBOL.y" ## Rename any existing SYMBOL column
  
  mapped_ids <- merge(data.frame(features[,unique_id_col]), mapping_table,by.x=1, by.y=1, all.x=T)
  colnames(mapped_ids)[1] <- "UNIQUE.ID"
  
  if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
    
    ## If a key maps to multiple identifiers...
    
    if(keep_all){
      
      ## ...keep all identifiers, and collapse rows with one-to-many mappings:
      
      mapped_ids <- mapped_ids %>%
        dplyr::group_by(UNIQUE.ID) %>%
        dplyr::summarise(
          TEMP=paste(!!as.name(map_to), sep=" | ")
        ) %>%
        as.data.frame()
      
      colnames(mapped_ids)[grep("TEMP", colnames(mapped_ids))] <- map_to
      
    } else {

      ## ...select only one:
      
      mapped_ids <- mapped_ids %>% 
        dplyr::group_by(UNIQUE.ID) %>% 
        dplyr::slice(1) %>%
        as.data.frame()
      
    }
    
  } ## if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
  
  ## Enforce that data is returned in its original order:
  
  mapped_ids <- mapped_ids[match(features[,unique_id_col], mapped_ids$UNIQUE.ID),]
  
  if(!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)){
    stop("!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)")
  }
  
  features <- data.frame(UNIQUE.ID=features1, TEMP=mapped_ids[,2], features[,-c(unique_id_col)])
  
  colnames(features)[2] <- map_to
  
  return(features)
  
}

mapAlias2Symbol <- function(features, unique_id_col, tables_path="/home/rebecca/omicon/mapping_tables/2022-11-02", keep_all=F, fill_NAs=F){
  
  mapping_table <- fread(list.files(pattern="ALIAS", path=file.path(tables_path, "homo_sapiens_mapping_tables/SYMBOL"), full.names=T), data.table=F)
  mapping_table <- as.data.frame(tidyr::separate_rows(mapping_table, ALIAS, sep=" \\| "))
  mapping_table[,2] <- toupper(mapping_table[,2])
  
  features1 <- features[,unique_id_col]
  features[,unique_id_col] <- toupper(features[,unique_id_col])
  
  mapped_ids <- merge(data.frame(features[,unique_id_col]), mapping_table, by.x=1, by.y=2, all.x=T)
  colnames(mapped_ids)[1] <- "UNIQUE.ID"
  
  if(fill_NAs){
    mapped_ids$SYMBOL[is.na(mapped_ids$SYMBOL)] <- mapped_ids$UNIQUE.ID[is.na(mapped_ids$SYMBOL)]
  }

  if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
    
    ## If a key maps to multiple identifiers...
    
    if(keep_all){
      
      ## ...keep all identifiers, and collapse rows with one-to-many mappings:
      
      mapped_ids <- mapped_ids %>% 
        dplyr::group_by(UNIQUE.ID) %>%
        dplyr::summarise(
          SYMBOL=paste(SYMBOL, collapse=" | ")
        ) %>%
        as.data.frame()
      
    } else { 
      
      ## ...select only one:
      
      mapped_ids <- mapped_ids %>%
        dplyr::group_by(UNIQUE.ID) %>%
        dplyr::slice(1)
      
    } 
    
  } ## if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
  
  ## Enforce that data is returned in its original order:
  
  mapped_ids <- mapped_ids[match(features[,unique_id_col], mapped_ids$UNIQUE.ID),]
  
  if(!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)){
    stop("!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)")
  }
  
  features <- data.frame(UNIQUE.ID=features1, SYMBOL=mapped_ids$SYMBOL, features[,-c(unique_id_col)])
  
  return(features)
  
}
