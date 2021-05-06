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

### this value came from RF_1_300.R
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
tuneGrid <- expand.grid(.mtry = best_mtry)

##this value came from RF_2_300.R 
max_node = 50


#identify best num of trees
store_maxtrees <- list()
for (ntree in c(100, 350, 500, 650, 750, 800, 850, 1000, 1500)) {
    set.seed(100)
    rf_maxtrees <- train(phenotype~.,
        data = train,
        method = "rf",
        metric = "Accuracy",
        tuneGrid = tuneGrid,
        trControl = control,
        importance = TRUE,
        nodesize = 3,
        maxnodes = max_node,
        ntree = ntree)
    key <- toString(ntree)
    store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
#best_ntree <- 250
##uncomment when you have the best tree.

##save to rData
save(results_tree, max_node, best_mtry, control, file = "data_mtry_part3_800t300g.RData")
