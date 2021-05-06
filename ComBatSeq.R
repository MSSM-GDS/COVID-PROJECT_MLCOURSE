.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") )
.libPaths()

library("BiocManager")
#BiocManager::install("sva", lib = "/hpc/users/jayarp02/R/lib/4.0.3/")

libs<-c("WGCNA","sva","corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR","AnnotationDbi","RColorBrewer","sva","GO.db","fitdistrplus","ff","plyranges","annotables","Rsamtools","GenomicFeatures","ggpubr","pheatmap","rms","dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace","reshape2","rmarkdown","org.Hs.eg.db","treemapify", "GSEABase", "varhandle","variancePartition", "ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges","scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize","colorspace","Vennerable","enrichR","cowplot","data.table","GSVA","caret","pROC","ROCR","mlbench","caretEnsemble","gbm","rpart","UpSetR","AUCell","RcisTarget","plyr","tidyverse","hrbrthemes","fmsb","colormap","viridis","survminer","survival","ggalluvial")



lapply(libs, require, character.only = TRUE)
rm(libs)

#actual datafile
data <- read.csv(file = "FINAL_INPUT_RAW.csv", header=TRUE, sep = ",")
#remove everything except gene counts
counts <- data[,]

print("counts")
ncol(counts)
print("genes column in dataframe")

#rename the row as gene names
rownames(counts) = data$genes

print("rownames of counts is the genes")
head(rownames(counts))

###remove GENES column - PJ
counts$genes <- NULL

sample_names = names(counts)

#add label
phenotype = read.csv("MASTER_PHENOTYPES_sorted.csv", header=TRUE, sep = ",")
y <- phenotype$Disease[match(colnames(counts), phenotype$GSEID_PatID)] #make sure the order is right
diseases = phenotype$Disease[match(colnames(counts), phenotype$GSEID_PatID)] #make sure the order is right
y <- as.numeric(as.factor(y))


###numeric vector equivalent for batch ID.
z <- phenotype$GSEID[match(colnames(counts), phenotype$GSEID_PatID)]
z <- as.numeric(as.factor(z))
batchid <- phenotype$GSEID[match(colnames(counts), phenotype$GSEID_PatID)]


sex <- phenotype$Sex[match(colnames(counts), phenotype$GSEID_PatID)]
s <- as.numeric(as.factor(sex))

condition <- phenotype$Condition[match(colnames(counts), phenotype$GSEID_PatID)]
c <- as.numeric(as.factor(condition))


pca_uncorrected_obj = prcomp(counts[,sample_names])
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)

pca_uncorrected[,"disease"] = diseases
pca_uncorrected[,"sex"] = sex
pca_uncorrected[,"batchid"] = batchid

pdf("before_ComBatSeq_pca.pdf", width = 10, height = 8)
###could do disease or condition or sex or group
ggplot(pca_uncorrected) + geom_point(aes(x=PC1, y=PC2, color = batchid, shape = disease))
dev.off()


corrected_data = ComBat_seq(counts = as.matrix(counts[,sample_names]), batch = z, group = y)
dim(counts)
dim(corrected_data)
write.csv(corrected_data,"x.BatchComBatSeq.csv")

###PCA corrected_data
pca_corrected_obj = prcomp(corrected_data[,sample_names])
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)
pca_corrected[,"disease"] = diseases
pca_corrected[,"sex"] = sex
pca_corrected[,"batchid"] = batchid

pdf("after_ComBatSeq_pca.pdf", width = 10, height = 8)
###could do disease or condition or sex or group
ggplot(pca_corrected) + geom_point(aes(x=PC1, y=PC2, color = batchid, shape = disease))
dev.off()
