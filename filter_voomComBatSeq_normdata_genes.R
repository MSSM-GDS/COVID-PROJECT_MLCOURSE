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

### read voom csv filte
voomComBatSeq_normcounts = read.csv(file = "x.BatchComBatSeq.csv", header=TRUE, sep = ",")
colnames(voomComBatSeq_normcounts)[1] <- "gene"
RFE_filteredgenes = read.table(file = "RFE_genes.txt", header=FALSE)
colnames(RFE_filteredgenes)[1] <- "gene"


filtered_data_csv = merge(voomComBatSeq_normcounts, RFE_filteredgenes, by.y="gene", by.x="gene", all.x = FALSE)
write.csv(filtered_data_csv, "VOOM_COMBATSEQ_NORM_FILTERED.csv")

filtered_data_csv = read.csv(file = "VOOM_COMBATSEQ_NORM_FILTERED.csv", header=TRUE, sep = ",")
filtered_data_csv$X = NULL
dim(filtered_data_csv)
nrow(filtered_data_csv)
write.csv(filtered_data_csv, "VOOM_COMBATSEQ_NORM_FILTERED.csv", row.names = FALSE)


### ensembl gene IDs to HGNC gene symbols
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mapped_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                                 filters = "ensembl_gene_id", values = as.vector(filtered_data_csv$gene),
                                 mart = mart))
#listAttributes(mart)
mapped_genes
write.csv(mapped_genes, "RFE_genes_symbols.csv")
