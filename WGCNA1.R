.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") )
.libPaths()

library("BiocManager")
libs<-c("WGCNA","sva","corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR","AnnotationDbi","RColorBrewer","sva","GO.db","fitdistrplus","ff","plyranges","annotables","Rsamtools","GenomicFeatures","ggpubr","pheatmap","rms","dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace","reshape2","rmarkdown","org.Hs.eg.db","treemapify", "GSEABase", "varhandle","variancePartition", "ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges","scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize","colorspace","Vennerable","enrichR","cowplot","data.table","GSVA","caret","pROC","ROCR","mlbench","caretEnsemble","gbm","rpart","UpSetR","AUCell","RcisTarget","plyr","tidyverse","hrbrthemes","fmsb","colormap","viridis","survminer","survival","ggalluvial")



lapply(libs, require, character.only = TRUE)
rm(libs)

options(stringsAsFactors = FALSE)
femData = read.csv(file = "x.BatchComBatSeq.csv", header=TRUE, sep = ",")
dim(femData)
head(femData[1:10, 1:10])
colnames(femData)[1] <- "gene"
datExpr0 = as.data.frame(t(femData[, -1])) * 1.0
names(datExpr0) = femData$gene
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "WGCNA_Clustering1.pdf", width = 400, height = 100)
par(cex = 0.9)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
clust = cutreeStatic(sampleTree, cutHeight = 10000000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
traitData = read.csv("MASTER_PHENOTYPES_sorted.csv", header=TRUE, sep = ",")
names(traitData)
head(traitData)
allTraits = traitData[, -1]
#allTraits = traitData[, -c(1:4, 7)]
dim(allTraits)
names(allTraits)
femaleSamples = rownames(datExpr)
traitRows = match(femaleSamples, allTraits$GSEID_PatID)
head(traitRows)
datTraits = allTraits[traitRows, -5]
rownames(datTraits) = allTraits[traitRows, 5]
collectGarbage();
sampleTree2 = hclust(dist(datExpr), method = "average")
datTraits2 <- data.matrix(datTraits)
traitColors = numbers2colors(datTraits2, signed = FALSE)
pdf(file = "WGCNA_Clustering_posttrim.pdf", width = 400, height = 100)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
dev.off()
save(datExpr, datTraits, datTraits2, file = "WGCNA-01-dataInput.RData")
