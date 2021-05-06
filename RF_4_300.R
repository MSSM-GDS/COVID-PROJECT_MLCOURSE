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


#save(data, train, test, rf_mtry, best_mtry, file = "data_mtry_part1_800t300g.RData")
#save(tuneGrid, results_mtry, control, file = "data_mtry_part2_800t300g.RData")
#save(results_tree, max_node, file = "data_mtry_part3_800t300g.RData")
load("data_mtry_part1_800t300g.RData")
load("data_mtry_part2_800t300g.RData")
load("data_mtry_part3_800t300g.RData")

### read voom csv file. Use Batch Normalized data. 
#voomed_counts = read.csv(file = "x.voom_afterBatchNorm.csv", header=TRUE, sep = ",")
#rownames(voomed_counts) <- voomed_counts$X
#voomed_counts$X <- NULL

##transpose
#data <- t(voomed_counts) #colnames should be genes and rownames should be samples; TRANSPOSED
#colnames(data) <- ev_markers$gene_symbol[match(ev_markers$entrez_id,colnames(data))]
#data <- data.frame(data)

#data <- as.data.frame(t(as.matrix(voomed_counts)))
#data2 = data
#dim(data)
#data = data[, 1:3000]

#add label
#data$phenotype <- phenotype$condition[match(rownames(data), phenotype$sample)] #make sure the order is right
#phenotype = read.csv("MASTER_PHENOTYPES_sorted.csv", header=TRUE, sep = ",")
#data$phenotype <- phenotype$Disease[match(rownames(data), phenotype$GSEID_PatID)] #make sure the order is right
#dim(data)

#best_mtry = 90

#i <- sample(nrow(data), 0.8*nrow(data), replace = FALSE)
#train <- data[i,]
#table(train$phenotype)
#bacterial     covid   healthy    others     viral
#      194       271       480       293       236
#test <- data[-i,]
#table(test$phenotype)
#bacterial     covid   healthy    others     viral
#       54        65       135        57        58

#control
#control <- trainControl(method="repeatedcv", number=10, repeats=10)

#tuneGrid <- expand.grid(.mtry = best_mtry)

#max_node = 75

##extracted from RF_3_300.R manually selected 500. 
#best_ntree = 500
best_ntree = 100

####final model.
fit_rf <- train(phenotype~.,
    data = train,
    method = "rf",
    metric = "Accuracy",
    tuneGrid = tuneGrid,
    trControl = control,
    importance = TRUE,
    nodesize = 5,
    ntree = best_ntree,
    maxnodes = max_node)

#evaluate model
prediction <- predict(fit_rf, test)
cm = confusionMatrix(prediction, as.factor(test$phenotype))
#Confusion Matrix and Statistics
#                     Reference
#Prediction            Primary Tumor Solid Tissue Normal
#  Primary Tumor                 101                   5
#  Solid Tissue Normal             2                   2
#               Accuracy : 0.9364
#                 95% CI : (0.8733, 0.974)
#    No Information Rate : 0.9364
#    P-Value [Acc > NIR] : 0.5988
#                  Kappa : 0.3328
# Mcnemar's Test P-Value : 0.4497
#            Sensitivity : 0.9806
#            Specificity : 0.2857
#         Pos Pred Value : 0.9528
#         Neg Pred Value : 0.5000
#             Prevalence : 0.9364
#         Detection Rate : 0.9182
#   Detection Prevalence : 0.9636
#      Balanced Accuracy : 0.6331
#
#       'Positive' Class : Primary Tumor

save(best_ntree, fit_rf, cm, file = "data_mtry_part4_800.RData")
