library(GSEABase)
library(parallel)
library(data.table)

source("/home/rebecca/code/map_identifiers/map_identifiers_function.R")

prep_MO_sets <- function(set_dir, MO_legend){

  set_files <- list.files(path=set_dir)
  set_files <- set_files[set_files%in%paste0(MO_legend$SetID, ".csv")]
  set_ids <- gsub(".csv", "", set_files)
  
  MO_legend <- MO_legend[match(set_ids, MO_legend$SetID),]
  
  if(!identical(set_ids, MO_legend$SetID)){
    stop("Some gene sets in the legend file were not found in the gene set directory, or vice versa.")
  }
  
  MO_sets <- lapply(1:length(set_files), function(i){
    read.csv(file.path(set_dir, set_files[i]), header=F)[,1]})
  
  if(n_distinct(MO_legend$SetID)!=nrow(MO_legend)){
    stop("Set identifiers are not unique!")
  }
  
  names(MO_sets) <- MO_legend$SetID
  
  save(MO_sets, MO_legend, file="/home/rebecca/gene_sets/MO/MO_sets.RData")
  
}

map_sets <- function(gene_sets, mapping_tables_dir, n_threads){
  
  ## For gene_sets with entries with multiple identfiers, select first value:
  gene_sets <- lapply(gene_sets, function(x) sapply(strsplit(x, "///", fixed=T), "[", 1))
  gene_sets <- lapply(gene_sets, function(x) sapply(strsplit(x, ";", fixed=T), "[", 1))
  ## Clean up:
  gene_sets <- lapply(gene_sets, function(x) x[x!="---"])
  gene_sets <- lapply(gene_sets, trimws)
  ## Map to most recent gene symbol:
  sets_mapped <- mclapply(gene_sets, function(x){
    mapAlias2Symbol(
      features=data.frame(x), 
      unique_id_col=1, 
      mapping_tables_dir, 
      keep_all_mappings=F, 
      fill_NAs=T
      )[,2]
    }, mc.cores=n_threads)
  return(sets_mapped)
  
}

get_broad_info <- function(idx, gene_sets){
  
  set <- gene_sets[[idx]]
  
  broad_info <- c(
    SetID=GSEABase::setIdentifier(set),
    SetName=GSEABase::setName(set),
    Category=GSEABase::collectionType(set)@category,
    SetSize=length(GSEABase::geneIds(set)),
    Species=GSEABase::organism(set),
    PubMed=paste(GSEABase::pubMedIds(set), collapse=" | "),
    Description=GSEABase::description(set)
  )
  return(broad_info)
  
}

prep_broad_sets <- function(xml_file, version){
  
  bsets <- GSEABase::getBroadSets(xml_file)
  
  broad_legend <- lapply(1:length(bsets), get_broad_info, bsets)
  broad_legend <- data.frame(
    do.call(rbind, broad_legend)
    )
  broad_sets <- lapply(bsets, GSEABase::geneIds)
  
  if(n_distinct(broad_legend$SetID)!=nrow(broad_legend)){
    stop("Set identifiers are not unique!")
  }
  
  names(broad_sets) <- broad_legend$SetID
  broad_sets <- broad_sets[broad_legend$SetSize>0]
  broad_legend <- broad_legend[broad_legend$SetSize>0,]
  broad_legend$SetSize <- as.numeric(broad_legend$SetSize)
  
  save(broad_sets, broad_legend, file=paste0(
    "/home/rebecca/gene_sets/broad/broad_sets_v", version, ".RData"))
  
}