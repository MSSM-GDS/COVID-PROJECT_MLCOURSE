##RFE
#x <- t(ev_voomed)
#phenotype$factor <- ifelse(phenotype$sample_type == "Primary_Tumor", 1, 0)
#y <- phenotype$factor[match(rownames(x), phenotype$sample)]
#y <- as.numeric(as.factor(y))
#identical(phenotype$sample, rownames(x)) #make sure sample order are the same
#number of attributes hoping to retain
#subsets <- seq(5, 10, 1)

.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") )
.libPaths()

library("BiocManager")
#‘ggradar’, ‘Vennerable’, ‘SCENIC’

#"papmap", "SCENIC",
libs<-c("corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR","AnnotationDbi","RColorBrewer","sva","GO.db","fitdistrplus","ff","plyranges","annotables","Rsamtools","GenomicFeatures","ggpubr","pheatmap","rms","dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace","reshape2","rmarkdown","org.Hs.eg.db","treemapify", "GSEABase", "varhandle","variancePartition", "ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges","scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize","colorspace","Vennerable","enrichR","cowplot","data.table","GSVA","caret","pROC","ROCR","mlbench","caretEnsemble","gbm","rpart","UpSetR","AUCell","RcisTarget","plyr","tidyverse","hrbrthemes","fmsb","colormap","viridis","survminer","survival","ggalluvial")

#install using BiocManager::install("") or install.packages("")
##for packages that arent available in R4.0.3 - use this method!
#packageURL <- "https://cran.r-project.org/src/contrib/Archive/tcR/tcR_2.3.2.tar.gz"
#install.packages(packageURL, repos=NULL, type="source", lib = "/hpc/users/jayarp02/R/lib/4.0.3/")


lapply(libs, require, character.only = TRUE)
rm(libs)

### read voom csv filte
voomed_counts = read.csv(file = "x.voom.csv", header=TRUE, sep = ",")
rownames(voomed_counts) <- voomed_counts$X
voomed_counts$X <- NULL

##transpose

x <- as.data.frame(t(as.matrix(voomed_counts)))
dim(x)

#add label
#data$phenotype <- phenotype$condition[match(rownames(data), phenotype$sample)] #make sure the order is right
phenotype = read.csv("MASTER_PHENOTYPES_sorted.csv", header=TRUE, sep = ",")

#phenotype$factor <- factor(ifelse(phenotype$sample_type == "Primary_Tumor", 1, 0))
#y <- phenotype$Disease[match(rownames(x), phenotype$sample)]
y <- phenotype$Disease[match(rownames(x), phenotype$GSEID_PatID)] #make sure the order is right
y <- as.numeric(as.factor(y))
#dim(y)

identical(phenotype$GSEID_PatID, rownames(x)) #make sure sample order are the same

#number of attributes hoping to retain
subsets <- seq(100, 800, 50)
#subsets <- seq(1, 10, 1)

set.seed(10)

#There are a number of pre-defined sets of functions for several models, including: linear regression (in the object lmFuncs), random forests (rfFuncs), naive Bayes (nbFuncs), bagged trees (treebagFuncs)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = TRUE,
                   number = 10)

print(ctrl)

Profile <- rfe(x, y,
                 sizes = subsets,
                 rfeControl = ctrl)

print(Profile)

pred = predictors(Profile)

save(Profile, pred, file = "Profile_RFE.RData")
