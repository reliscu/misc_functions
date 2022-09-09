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
  
  mapping_table <- fread(list.files(pattern=map_to, path=file.path(mapping_tables_dir, tables_dir), full.names=T), data.table=F)
  mapping_table[,1] <- toupper(mapping_table[,1])
  
  features_orig <- features[,unique_id_col]
  features[,unique_id_col] <- toupper(features[,unique_id_col])
  colnames(features)[colnames(features)=="SYMBOL"] <- "SYMBOL.y" ## Rename any existing SYMBOL column
  
  mapped_ids <- merge(features[,unique_id_col], mapping_table, by.x=1, by.y=1, all.x=T)
  colnames(mapped_ids)[1] <- "UNIQUE.ID"
  
  if(keep_all_mappings==F){
    
    ## If a key maps to multiple identifiers, select only one:
    
    if(length(duplicated(mapped_ids$UNIQUE.ID))>0){
      
      mapped_ids <- mapped_ids %>% 
        dplyr::group_by(UNIQUE.ID) %>% 
        dplyr::slice(1)
      
      # tidyr::separate_rows(!!as.name(map_to), sep=" \\| ") %>%
      #   dplyr::arrange(!!as.name(map_to)) %>%
      
    }
    
  } else {
    
    # mapped_ids <- mapped_ids %>% 
    #   dplyr::group_by(UNIQUE.ID) %>% 
    #   dplyr::summarise(TEMP=paste(!!as.name(map_to), sep=" | ")) %>%
    #   as.data.frame()
    # colnames(mapped_ids)[grep("TEMP", colnames(mapped_ids))] <- map_to
    
  }
  
  ## Order of identifiers must match input order:
  
  mapped_ids <- mapped_ids[match(features[,unique_id_col], mapped_ids$UNIQUE.ID),]
  
  if(!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)){
    stop("!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)")
  }
  
  features <- data.frame(UNIQUE.ID=features_orig, TEMP=mapped_ids[,2], features[,-c(unique_id_col)])
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
  
  mapping_table <- fread(list.files(pattern="ALIAS", path=file.path(mapping_tables_dir, hs_dir, "SYMBOL"), full.names=T), data.table=F)
  mapping_table <- data.frame(separate_rows(mapping_table, ALIAS, sep=" \\| "))
  mapping_table[,2] <- toupper(mapping_table[,2])
  
  features_orig <- features[,unique_id_col]
  features[,unique_id_col] <- toupper(features[,unique_id_col])
  
  mapped_ids <- merge(
    data.frame(features[,unique_id_col]), 
    mapping_table, by.x=1, by.y=2, all.x=T
    )
  
  colnames(mapped_ids)[1] <- "UNIQUE.ID"
  
  if(fill_NAs==T){
    mapped_ids$SYMBOL[is.na(mapped_ids$SYMBOL)] <- mapped_ids$UNIQUE.ID[is.na(mapped_ids$SYMBOL)]
  }
  
  ## For aliases that map to 2+ symbols:
  
  if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
    
    if(keep_all_mappings==F){
      
      mapped_ids <- mapped_ids %>%
        dplyr::group_by(UNIQUE.ID) %>%
        dplyr::slice(1)
      
    } else { ## if(keep_all_mappings==F){
      
      mapped_ids <- mapped_ids %>% 
        dplyr::group_by(UNIQUE.ID) %>%
        dplyr::summarise(SYMBOL=paste(SYMBOL, collapse=" | "))
      
    } ## if(keep_all_mappings==F){} else {
    
  } ## if(nrow(dupl_ids)>0){
  
  mapped_ids <- mapped_ids[match(features[,unique_id_col], mapped_ids$UNIQUE.ID),]
  
  if(!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)){
    stop("!identical(features[,unique_id_col], mapped_ids$UNIQUE.ID)")
  }
  
  features <- data.frame(
    UNIQUE.ID=features_orig, 
    SYMBOL=mapped_ids$SYMBOL, 
    features[,-c(unique_id_col)]
  )
  
  return(features)
  
}


#dupl_ids <- mapped_ids[mapped_ids$UNIQUE.ID%in%mapped_ids$UNIQUE.ID[duplicated(mapped_ids$UNIQUE.ID)],]

# if(sum(duplicated(mapped_ids$UNIQUE.ID))>0){
#   
#   mapped_ids <- mapped_ids %>%
#     dplyr::group_by(UNIQUE.ID) %>%
#     dplyr::slice(1)
#   
#   # nondupl_ids <- mapped_ids[!mapped_ids$UNIQUE.ID%in%dupl_ids$UNIQUE.ID,]
#   # 
#   # ## If ALIAS==SYMBOL, choose that symbol:
#   # 
#   # alias_symbol_match <- dupl_ids[dupl_ids$UNIQUE.ID==dupl_ids$SYMBOL,]
#   # alias_symbol_mismatch <- dupl_ids[!dupl_ids$UNIQUE.ID%in%alias_symbol_match$UNIQUE.ID,]
#   # 
#   # if(nrow(alias_symbol_mismatch)>0){
#   #   
#   #   ## For remaining aliases, keep the first symbol:
#   #   
#   #   alias_symbol_mismatch <- alias_symbol_mismatch %>%
#   #     dplyr::arrange(SYMBOL) %>%
#   #     dplyr::group_by(UNIQUE.ID) %>%
#   #     dplyr::slice(1)
#   
# } ## if(nrow(dupl_ids)>0){
# 
# #mapped_ids <- rbind(nondupl_ids, alias_symbol_match, alias_symbol_mismatch)
