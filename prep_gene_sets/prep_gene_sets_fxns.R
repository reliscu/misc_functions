library(GSEABase)
library(data.table)
library(future.apply)
library(openblasctl)

openblas_set_num_threads(10)
options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

source("/home/rebecca/code/misc/map_identifiers/map_identifiers_function.R")

prep_MO_sets <- function(projectname, set_dir, MO_legend, out_dir){

  set_files <- list.files(path=set_dir)
  set_files <- set_files[is.element(set_files, paste0(MO_legend$SetID, ".csv"))]
  set_ids <- gsub(".csv", "", set_files)
  
  MO_legend <- MO_legend[match(set_ids, MO_legend$SetID),]
  
  if(!identical(set_ids, MO_legend$SetID)){
    stop("Some gene sets in the legend file were not found in the gene set directory, or vice versa.")
  }
  
  MO_sets <- lapply(1:length(set_files), function(i){
    read.csv(file.path(set_dir, set_files[i]), header=F)[,1]
  })

  if(n_distinct(MO_legend$SetID)!=nrow(MO_legend)){
    stop("Set identifiers are not unique!")
  }
  
  names(MO_sets) <- MO_legend$SetID
  
  save(MO_sets, MO_legend, file=paste0(out_dir, "/", projectname, "_sets.RData"))
  
}

map_sets <- function(projectname, gene_sets, legend, mapping_tables_dir, out_dir, n_threads){
  
  ## For gene sets with entries with multiple entires per identifier, select first entry:
  
  gene_sets <- lapply(gene_sets, function(x) sapply(strsplit(x, "///", fixed=T), "[", 1))
  gene_sets <- lapply(gene_sets, function(x) sapply(strsplit(x, ";", fixed=T), "[", 1))
  
  ## Clean up:
  
  gene_sets <- lapply(gene_sets, function(x) x[x!="---"])
  gene_sets <- lapply(gene_sets, trimws)
  
  ## Map to most recent gene symbol:
  
  gene_sets_mapped <- future_lapply(gene_sets, function(set){
    
    return(
      mapAlias2Symbol(
        features=data.frame(set), 
        unique_id_col=1, 
        mapping_tables_dir, 
        keep_all_mappings=F, 
        fill_NAs=T
      )[,2]
    )
    
  })
  
  if(!identical(names(gene_sets), legend$SetID)){
    stop("!identical(names(gene_sets), legend$SetID)")
  }
  
  names(gene_sets_mapped) <- legend$SetID
  
  save(gene_sets_mapped, legend, file=paste0(out_dir, "/", projectname, "_sets_mapped.RData"))
  
}

get_broad_info <- function(set){

  return(c(
    SetID=GSEABase::setIdentifier(set),
    SetName=GSEABase::setName(set),
    Category=GSEABase::collectionType(set)@category,
    SetSize=length(GSEABase::geneIds(set)),
    Species=GSEABase::organism(set),
    PubMed=paste(GSEABase::pubMedIds(set), collapse=" | "),
    Description=GSEABase::description(set)
  ))

}

prep_broad_sets <- function(projectname, xml_file, version, out_dir){
  
  bsets <- GSEABase::getBroadSets(xml_file)
  
  broad_legend <- lapply(bsets, get_broad_info)
  broad_legend <- data.frame(do.call(rbind, broad_legend))
  broad_sets <- lapply(bsets, GSEABase::geneIds)
  
  if(n_distinct(broad_legend$SetID)!=nrow(broad_legend)){
    stop("Set identifiers are not unique!")
  }
  
  broad_sets <- broad_sets[broad_legend$SetSize>0]
  broad_legend <- broad_legend[broad_legend$SetSize>0,]
  broad_legend$SetSize <- as.numeric(broad_legend$SetSize)
  
  names(broad_sets) <- broad_legend$SetID
  
  save(broad_sets, broad_legend, file=paste0(out_dir, "/", projectname, "_", version, ".RData"))
  
}