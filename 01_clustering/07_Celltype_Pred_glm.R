# ====================================================================
# Author: Zhou Zhe
# Function: Use lasso to predict celltype based on metabolic genes expression
# Groups: TpN means using all cells (tumor + normal), Tumor means using only tumor cells
#         resample_cells group means subsample cells directly, the other group with no prefix
#         means dividing samples into train and test first, and then subsample cells
#         *Finally, use no prefix TpN group, which train dataset and test dataset including 
#         different samples*
# Version: 1.0
# Date: Feb 21, 2023
# ====================================================================


rm(list = ls())
library(Seurat)
library(dplyr)
library(glmnet)
library(ROCR)
require(doMC)
library(plyr)
library(ggplot2)
library(scMetab)
library(patchwork)

registerDoMC(cores = 12)
seed.use = 2021
dir.used <- "TpN/"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))

# I. GLM ===========================================================================================
sink(paste0("data/rds/Celltype_Prediction/", dir.used, "TpN.log"))
for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste(Sys.time(), tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(paste0("data/rds/MetabGene/", tumor, ".TpN.unintegrated.rds"))
    # seu_obj <- subset(seu_obj, TorN == "T")
    print(table(seu_obj$celltype))
    
    ## features and targets
    x <- t(GetAssayData(seu_obj, slot = "scale.data"))
    used_celltypes <- sort(unique(seu_obj$celltype))
    
    # subsample samples and cells
    nCells <- table(seu_obj$patientID)
    set.seed(seed.use)
    trainSam <- sample(names(nCells), (length(nCells)+1)/2)
    testSam <- setdiff(names(nCells), trainSam)
    print(paste0(sum(nCells[trainSam]), " cells in train sample and ", sum(nCells[testSam]), " cells in test sample!"))
    trainCells <- sample(colnames(seu_obj)[seu_obj$patientID %in% trainSam], min(10000, sum(nCells[trainSam])))
    testCells <- sample(colnames(seu_obj)[seu_obj$patientID %in% testSam], min(10000, sum(nCells[testSam])))
    set.seed(NULL)
    
    # set.seed(seed.use)
    # subCells <- sample(colnames(seu_obj), min(20000, dim(seu_obj)[2]))
    # mid <- ceiling(length(subCells)/2)
    # trainCells <- subCells[1:mid]
    # testCells <- subCells[(mid+1):length(subCells)]
    # set.seed(NULL)
    
    roc.list = list()
    tmpfile <- paste0("plot/1_clustering/Celltype_Pred/", dir.used, tumor, "_ROC.pdf")
    pdf(tmpfile, width = 7, height = 7)
    for (tarCT in used_celltypes) {
        
        cat(paste(Sys.time(), "--------", tarCT, "--------\n"))
        y <- ifelse(seu_obj$celltype == tarCT, TRUE, FALSE)
        
        ## cross validation
        system.time(cvfit <- cv.glmnet(x[trainCells, ], y[trainCells], family = "binomial", type.measure = "auc", parallel = T))
        plot(cvfit, main = paste(tarCT, "para-tuning CV"))
        # c(cvfit$lambda.min, cvfit$lambda.1se)
        
        # coefficient matrix
        coef1 <- coef(cvfit$glmnet.fit, s = cvfit$lambda.min, exact = F)
        coef2 <- coef(cvfit$glmnet.fit, s = cvfit$lambda.1se, exact = F)
        coef.df <- data.frame(cbind(coef1, coef2))
        
        ## confusion matrix of train and test data
        a <- predict(cvfit$glmnet.fit, s = cvfit$lambda.min, newx = x[trainCells, ], type = "class")
        print(table(a, y[trainCells]))
        b <- predict(cvfit$glmnet.fit, s = cvfit$lambda.min, newx = x[testCells, ], type = "class")
        print(table(b, y[testCells]))
        
        ## ROC and precision-recall curve
        prob.train <- predict(cvfit$glmnet.fit, s = cvfit$lambda.min, newx = x[trainCells, ], type = "response")
        pred.train <- prediction(prob.train, y[trainCells])
        roc.train <- performance(pred.train, "tpr", "fpr")
        plot(roc.train, main = paste(tarCT, "training data"))
        auc.train <- performance(pred.train, "auc")@y.values[[1]]
        
        prob.test <- predict(cvfit$glmnet.fit, s = cvfit$lambda.min, newx = x[testCells, ], type = "response")
        pred.test <- prediction(prob.test, y[testCells])
        roc.test <- performance(pred.test, "tpr", "fpr")
        plot(roc.test, main = paste(tarCT, "test data"))
        auc.test <- performance(pred.test, "auc")@y.values[[1]]
        
        roc.list[[tarCT]] <- list(roc.train = roc.train, auc.train = auc.train, 
                                  roc.test = roc.test, auc.test = auc.test)
        # prc <- performance(pred, "prec", "rec")
        # plot(prc)
    }
    dev.off()
    
    # save roc.list
    saveRDS(roc.list, file = paste0("data/rds/Celltype_Prediction/", dir.used, tumor, ".pred.rds"))
}
sink()


# II. Plot =========================================================================================


pdf(paste0("plot/1_clustering/Celltype_Pred/", dir.used, "Pooled_ROC.pdf"), width = 3.3, height = 1.8)
tmp.text <- element_text(family="sans", size=6)
for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    roc <- readRDS(paste0("data/rds/Celltype_Prediction/", dir.used, tumor, ".pred.rds"))
    mycolor <- celltype_color[names(celltype_color) %in% names(roc)]
    roc <- roc[names(mycolor)]
    
    roc.train <- ldply(names(roc), function(celltype) {
        data.frame(fpr = roc[[celltype]]$roc.train@x.values[[1]], 
                   tpr = roc[[celltype]]$roc.train@y.values[[1]], 
                   celltype = celltype)
    })
    roc.train$celltype <- factor(roc.train$celltype, levels = names(mycolor))
    auc.train <- sapply(roc, function(x) x$auc.train)
    p1 = ggplot(roc.train, aes(x = fpr, y = tpr, group = celltype)) + 
        geom_line(aes(color = celltype), linewidth = 0.25) + 
        scale_color_manual(values = mycolor, labels = paste(names(auc.train), sprintf("%0.3f", auc.train))) + 
        ggtitle("Training data") + xlab("1 - specificity") + ylab("Sensitivity") + 
        theme_bw()+
        theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text,
              plot.title = element_text(hjust = 0.5, size = 6), 
              panel.grid = element_blank(), legend.position = c(0.6,0.5), 
              legend.key.size = unit(1, "mm"), legend.title = element_blank())
    
    roc.test <- ldply(names(roc), function(celltype) {
        data.frame(fpr = roc[[celltype]]$roc.test@x.values[[1]], 
                   tpr = roc[[celltype]]$roc.test@y.values[[1]], 
                   celltype = celltype)
    })
    roc.test$celltype <- factor(roc.test$celltype, levels = names(mycolor))
    auc.test <- sapply(roc, function(x) x$auc.test)
    p2 = ggplot(roc.test, aes(x = fpr, y = tpr, group = celltype)) + 
        geom_line(aes(color = celltype), linewidth = 0.25) + 
        scale_color_manual(values = mycolor, labels = paste(names(auc.test), sprintf("%0.3f", auc.test))) + 
        ggtitle("Test data") + xlab("1 - specificity") + ylab("Sensitivity") + 
        theme_bw()+
        theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text,
              plot.title = element_text(hjust = 0.5, size = 6), 
              panel.grid = element_blank(), legend.position = c(0.6,0.5), 
              legend.key.size = unit(1, "mm"), legend.title = element_blank())
    
    print(p1 + p2 + plot_annotation(tumor, theme = theme(plot.title = element_text(size = 6))))
}
dev.off()
