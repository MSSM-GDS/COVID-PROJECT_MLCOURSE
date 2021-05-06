.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") )
.libPaths()

#library("BiocManager")
#BiocManager::install(version = "3.12", repo = "http://cran.us.r-project.org", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
#remotes::install_version("RSQLite", version = "2.2.5", repo = "http://cran.us.r-project.org", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
#BiocManager::install(c("Biobase", "GenomicRanges","org.Hs.eg.db","GO.db","GSEABase","limma","edgeR","AnnotationDbi","enrichR","GSVA", "DESeq2","Rfast","annotables", "snm", "doMC", "tibble", "gbm"), repo = "http://cran.us.r-project.org", lib="/hpc/users/jayarp02/R/lib/4.0.3/")

#install.packages("tidyverse", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
#install.packages("ggplot2", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
#install.packages("ggrepel", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
#install.packages("RColorBrewer", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
#install.packages(c("ggpubr","ggbeeswarm","tibble","rmarkdown","pheatmap","colorspace","colormap","reshape2","plyr","data.table","rpart","scales","fitdistrplus","plyranges","viridis","doParallel","factoextra"), lib="/hpc/users/jayarp02/R/lib/4.0.3/")


#################
## STAGE 1: DATA QC AND NORMALIZATION
#################

libs<-c("tidyr","dplyr","ggplot2","ggrepel","RColorBrewer","ggpubr","ggbeeswarm",
        "tibble","rmarkdown","pheatmap","colorspace","colormap","reshape2","plyr",
        "data.table","rpart","scales","fitdistrplus","plyranges","tidyverse","viridis",
        "doParallel","factoextra","org.Hs.eg.db","GO.db","GSEABase",
        "limma","edgeR","AnnotationDbi","enrichR","GSVA", "DESeq2",
        "Rfast","annotables", "snm", "doMC", "tibble", "gbm")
lapply(libs, require, character.only = TRUE)
rm(libs)

################
##THINGS TO IMPORT PRIOR TO RUNNING
################
#phenotype file
phenotype <- read.csv(file = "MASTER_PHENOTYPES_sorted.csv", header=TRUE, sep = ",")
#change according to phenotype file label

group <- as.factor(phenotype$Condition)
print("group")
length(group)
#head(group)

#actual datafile
#data <- read.csv(file = "FINAL_INPUT_RAW.csv", header=TRUE, sep = ",")
data <- read.csv(file = "x.BatchComBatSeq.csv", header=TRUE, sep = ",")
#remove everything except gene counts
counts <- data[,]

print("counts")
ncol(counts)
#head(counts)

print("genes column in dataframe")
head(data$genes)

#rename the row as gene names
rownames(counts) = data$genes

print("rownames of counts is the genes")
head(rownames(counts))


###remove GENES column - PJ
counts$genes <- NULL
print("counts data frame without genes")
#head(counts)
ncol(counts)

samples <- data.frame(colnames(counts))
print("samples")
head(samples)


#####filter out zero expressed genes; determine by function
#keep.exprs <- filterByExpr(, group=group)
#x <- counts[keep.exprs,, keep.lib.sizes=FALSE]

####determine filtering manually; change 500 to number of samples in total
#table(rowSums(counts>30)>300)
#table(rowSums(counts>50)>300)
#table(rowSums(counts>50)>100)
table(rowSums(counts>0)>0)

keep.expr <- rowSums(counts>0)>0
x <- counts[keep.expr,]
print("dimensions of X:")
dim(x)
dcpmx <- DGEList(counts=x,group=group, genes = data[,1])

print("dcpmx done:")
#head(dcpmx)

################
##Normalization
################
#####log2 transformation
#x.log <- log2(x+1)
#x.log[1:5,1:5]
#write.csv(x.log, "x.log.csv")


#####TMM
#TMM <- calcNormFactors(dcpmx, method = "TMM")
#x.TMM <- cpm(TMM)
#x.TMM[1:5, 1:5]
#write.csv(x.TMM, "x.TMM.csv")

#####Voom
#create design matrix
designV <- model.matrix(~0 + Disease + GSEID + Sex , data = phenotype)
#colnames(designV) <- unique(group)
voom <- voom(dcpmx, design = designV, plot = TRUE, normalize.method="quantile")
x.voom <- voom$E
x.voom[1:5,1:5]
write.csv(x.voom, "x.voom_afterBatchNorm.csv")

#Voom + SNM
#bio.var <- model.matrix(~group + disease) #disease
#colnames(bio.var) <- unique(group)
#adj.var <- model.matrix(~ age + sex) #age/sex
#colnames(adj.var) <- unique(group)
#snm <- snm(raw.dat = x.voom,
#			bio.var = bio.var,
#			adj.var = adj.var,
#			rm.adj=TRUE,
#			verbose = TRUE,
#			diagnose = FALSE)
#x.snm <- t(snm$norm.dat)
#x.snm[1:5,1:5]
#write.csv(x.snm, "x.snm.csv")

#DESeq2
print("phenotype or colData")
#print(table(phenotype))

print("group or design")
print(table(group))

dds <- DESeqDataSetFromMatrix(countData = x, colData = phenotype, design = ~Condition + Sex) #binary labels
#print("dds")
print(head(dds))

dds <- estimateSizeFactors(dds)
x.deseq2 <- counts(dds, normalized=TRUE)
print(x.deseq2[1:5,1:5])

write.csv(x.deseq2,"x.deseq2.csv")
print("deseq2 output")

################
##PCA
################
#before
pca.count <- counts
pca <- prcomp(t(pca.count))
pca.file <- cbind(phenotype, pca$x)
pdf("before_pca_plot_CondSex_PC23.pdf", width = 10, height = 8)
###could do disease or condition or sex or group
ggplot(pca.file) + geom_point(aes(x=PC2, y=PC1, color = Condition, shape = Sex))
#+ ylim(NA, 5000000)
#+  xlim(NA, 1000000)
dev.off()

#log
#pca.count <- x.log
#pca <- prcomp(t(pca.count))
#pca.file <- cbind(phnotype, pca$pca.count)
#pdf("log_pca.plot", width = 10, height = 8)
#ggplot(pca.file) + geom_point(aes(x=PC1, y=PC2, color = condition,
#	shape = disease)) #could do disease as well
#dev.off()

#TMM
#pca.count <- x.TMM
#pca <- prcomp(t(pca.count))
#pca.file <- cbind(phnotype, pca$pca.count)
#pdf("TMM_pca.plot", width = 10, height = 10)
#ggplot(pca.file) + geom_point(aes(x=PC1, y=PC2, color = condition,
#	shape = disease)) #could do disease as well
#dev.off()

#voom
pca.count <- x.voom
pca <- prcomp(t(pca.count))
pca.file <- cbind(phenotype, pca$x)
pdf("voom_pca_postBatchNorm_plot_SexDiseaseGSE_PC13.pdf", width = 10, height = 8)
ggplot(pca.file) + geom_point(aes(x=PC1, y=PC3, color = GSEID, shape = Disease)) #could do disease as well
dev.off()


#snm
#pca.count <- x.snm
#pca <- prcomp(t(pca.count))
#pca.file <- cbind(phnotype, pca$pca.count)
#pdf("snm_pca.plot", width = 10, height = 8)
#ggplot(pca.file) + geom_point(aes(x=PC1, y=PC2, color = condition,
#	shape = disease)) #could do disease as well
#dev.off()

#DESeq2
pca.count <- x.deseq2
pca <- prcomp(t(pca.count))
pca.file <- cbind(phenotype, pca$x)
pdf("deseq_pca.plot_ConditionSex_PC12.pdf", width = 10, height = 10)
ggplot(pca.file) + geom_point(aes(x=PC2, y=PC1, color = Condition, shape = Sex)) #could do disease as well
dev.off()
