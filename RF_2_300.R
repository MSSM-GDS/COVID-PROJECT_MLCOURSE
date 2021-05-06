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

load("data_mtry_part1_800t300g.RData")
dim(data)

print(best_mtry)


####acquired this from outputs from 1
best_mtry = 95

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
control <- trainControl(method="repeatedcv", number=5, repeats=10)

####identify best maxnodes; tree depth - manual
#kappa = compares observed accuracy with expected accuracy;
	#not only evaluate a single classifier but amongst classifier themselves
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in seq(from=35, to=125, by=10)) {
    set.seed(100)
    rf_maxnode <- train(phenotype~.,
        data = train,
        method = "rf",
        metric = "Accuracy",
        tuneGrid = tuneGrid,
        trControl = control,
        importance = TRUE,
        nodesize = 70,
        maxnodes = maxnodes,
        ntree = 300)
    current_iteration <- toString(maxnodes)
    store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)
###max_node <- 5 - needs to be chosen manually. -PJ choose max or highest average. kappa values shoud not be low(0s or negatives)
###uncomment when you get max_node and set max_node.

save(tuneGrid, results_mtry, max_node, control, file = "data_mtry_part2_800t300g.RData")
