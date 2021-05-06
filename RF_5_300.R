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
#save(best_ntree, fit_rf, cm, file = "data_mtry_part4_800.RData")
load("data_mtry_part1_800t300g.RData")
load("data_mtry_part2_800t300g.RData")
load("data_mtry_part3_800t300g.RData")
load("data_mtry_part4_800t300g.RData")

dim(data)

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

####identify best maxnodes; tree depth - manual
#kappa = compares observed accuracy with expected accuracy;
	#not only evaluate a single classifier but amongst classifier themselves

#tuneGrid <- expand.grid(.mtry = best_mtry)

#max_node = 75

#best_ntree = 300

print(cm)


print(fit_rf)

###this will give actual feaures that you will use.
feature_varimp = varImp(fit_rf)
#rf variable importance
#        Importance
#NKX3.1     100.000
#SMYD3       99.733
#TET3        97.114
#CYLD        73.502
#SNORA54     66.577
#HAS2        48.142
#MXD4        24.572
#UPK1B       21.923
#FOXO3       19.612
#ESR1        19.475
#IRF1         3.842
#BRCA1        0.000


feature_varimp


#ROC -- change variables
library(pROC)
roc.predict <- predict(fit_rf, test, type = "prob")
roc.predict

#result.roc <- roc(test$phenotype, roc.predict$Disease, plot = T)
#result.roc

#pdf("FINAL_ROC_AUC_curve_300.pdf", height = 10, width = 8)
#plot(result.roc, print.auc=T)
##Area under the curve: 0.846
#dev.off()
