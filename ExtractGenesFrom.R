.libPaths( c( .libPaths(), "/Users/jayaramanp/workspace/R/lib/4.04/") )
.libPaths()

library("gdata")
library("tidyverse")
library("AnnotationDbi")
library("org.Hs.eg.db")
#library("ensembldb")
#library("EnsDb.Hsapiens.v86")
library("biomaRt")
library("data.table")
library("dplyr")
options(stringsAsFactors = FALSE)
rf4_fs_db <- read.csv(file = "RF4fs_varImpgenes.csv", header=TRUE, sep = ",")
rf4_fs_filtereddb = rf4_fs_db %>% filter_all(all_vars(. > 30))


### ensembl gene IDs to HGNC gene symbols
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mapped_genes_rf4fs <- data.frame(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                                 filters = "ensembl_gene_id", values = as.vector(rf4_fs_filtereddb$X),
                                 mart = mart))
#listAttributes(mart)
mapped_genes_rf4fs
write.csv(mapped_genes_rf4fs, "RF4fs_filteredgenesymbols.csv")

rfe_fs_db <- read.csv(file = "RFE_genes_symbols.csv", header=TRUE, sep = ",")
wgcna_genes_db <- read.csv(file = "WGCNA_genes.csv", header=FALSE, sep = ",")

list1 = as.vector(mapped_genes_rf4fs$ensembl_gene_id)
list2 = as.vector(rfe_fs_db$ensembl_gene_id)
list3 = as.vector(wgcna_genes_db$V1)
intersect = Reduce(intersect, list(list1, list2),list3)
intersect

### ensembl gene IDs to HGNC gene symbols
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mapped_genes_intersect <- data.frame(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                                       filters = "ensembl_gene_id", values = intersect,
                                       mart = mart))
#listAttributes(mart)
mapped_genes_intersect
write.csv(mapped_genes_intersect, "VennDiagram_geneSetIntersect.csv")
