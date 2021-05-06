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
voomed_counts = read.csv(file = "VOOM_COMBATSEQ_NORM_FILTERED.csv", header=TRUE, sep = ",")
rownames(voomed_counts) <- voomed_counts$gene
voomed_counts$gene <- NULL

##transpose
#data <- t(voomed_counts) #colnames should be genes and rownames should be samples; TRANSPOSED
#colnames(data) <- ev_markers$gene_symbol[match(ev_markers$entrez_id,colnames(data))]
#data <- data.frame(data)

data <- as.data.frame(t(as.matrix(voomed_counts)))
#data2 = data
dim(data)
#data = data[, 1:3000]

#add label
#data$phenotype <- phenotype$condition[match(rownames(data), phenotype$sample)] #make sure the order is right
phenotype = read.csv("MASTER_PHENOTYPES_sorted.csv", header=TRUE, sep = ",")
data$phenotype <- phenotype$Disease[match(rownames(data), phenotype$GSEID_PatID)] #make sure the order is right

set.seed(100)

##split into training and test dataset.
#set proportion
i <- sample(nrow(data), 0.8*nrow(data), replace = FALSE)
train <- data[i,]
table(train$phenotype)
#bacterial     covid   healthy    others     viral
#      194       271       480       293       236
test <- data[-i,]
table(test$phenotype)
#bacterial     covid   healthy    others     viral
#       54        65       135        57        58



#control
control <- trainControl(method="repeatedcv", number=5, repeats=5, search='grid', allowParallel = TRUE)
#training
set.seed(100)

#mtry <- sqrt(ncol(train)) #number of variables randomly sampled as attribute at each split
#tunegrid <- expand.grid(.mtry=mtry)
rf_untuned <- train(phenotype~.,
	data=train,
	method="rf",
	metric="Accuracy",
	trControl=control)
print(rf_untuned)
#Random Forest

#100 samples
# 15 predictor
#  2 classes: 'Primary Tumor', 'Solid Tissue Normal'

#No pre-processing
#Resampling: Cross-Validated (10 fold, repeated 10 times)
#Summary of sample sizes: 395, 397, 396, 396, 395, 395, ...
#Resampling results across tuning parameters:

#  mtry  Accuracy  Kappa
#   2    0.832     0.664
#   8    0.839     0.678
#  15    0.827     0.654

#Accuracy was used to select the optimal model using the largest value.
#The final value used for the model was mtry = 7.

#identify best mtry
set.seed(100)
tuneGrid <- expand.grid(.mtry = c(10, 25, 50, 75, 90, 95, 105,  115, 125))
rf_mtry <- train(phenotype~.,
    data = train,
    method = "rf",
    metric = "Accuracy",
    tuneGrid = tuneGrid,
    trControl = control,
    importance = TRUE,
    nodesize = 14,
    ntree = 500)
print(rf_mtry)
max(rf_mtry$results$Accuracy)
#store best mtry
best_mtry <- rf_mtry$bestTune$mtry
save(data, train, test, rf_mtry, best_mtry, file = "data_mtry_part1_800t300g.RData")
