source("/home/rebecca/code/prep_gene_sets/prep_sets.R")

## Format gene sets

xml_file <- "/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml"
version <- 7
out_dir <- "/home/rebecca/gene_sets/broad/"

prep_broad_sets(xml_file, version, out_dir)

set_dir <- "/home/shared/genesets/OurSets/"
oldham_legend <- read.csv("/home/shared/genesets/OurSets/oldham_legend.csv")
out_dir <- "/home/rebecca/gene_sets/oldham/"

prep_oldham_sets(set_dir, oldham_legend, out_dir)

## Map to latest identifiers

mapping_tables_dir <- "/home/rebecca/omicon/mapping_tables/"

load("/home/rebecca/gene_sets/oldham/oldham_sets.RData")
out_dir <- "/home/rebecca/gene_sets/oldham/"
sets <- oldham_sets
legend <- oldham_legend

oldham_sets_mapped <- map_sets(sets, legend, mapping_tables_dir, out_dir)

save(oldham_sets_mapped, oldham_legend, file=file.path(out_dir, "oldham_sets_mapped.RData"))

load("/home/rebecca/gene_sets/broad/broad_sets_v7.RData")
out_dir <- "/home/rebecca/gene_sets/broad/"
sets <- broad_sets
legend <- broad_legend

broad_sets_mapped <- map_sets(sets, legend, mapping_tables_dir, out_dir)

save(broad_sets_mapped, broad_sets, file=file.path(out_dir, "broad_sets_v7_mapped.RData"))