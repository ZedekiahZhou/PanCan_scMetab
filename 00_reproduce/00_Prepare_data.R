# ==================================================================================================
# Author: Zhou Zhe
# Program: Prepare seurat object, count cells for each sample and each cell type
# Version: 1.0
# Date: Jul 21, 2023
# ==================================================================================================

rm(list=ls())
library(Seurat)
library(scMetab)
library(fastSave)
# cl = makeCluster(8)


dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
dir.used <- "00_reproduce/Cell_counts_celltype_x_sample/"
rownames(dataset) <- dataset$DataSets

# Prepare metabolic gene expression matrix for NMF ========
for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(Sys.time(), " ------ ", tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(paste0("../", dataset$Directory[i], "/data/rds/", dataset$Prefix[i], ".TpN.merged.rds"))
    table(seu_obj$samID)
    
    if (!("samID" %in% colnames(seu_obj[[]]))) seu_obj$samID = seu_obj$patientID
    if (tumor == "PRAD") seu_obj$samID <- seu_obj$samID2
    
    table(seu_obj$sub_celltype)
    
    ## distinguish CD4+ and CD8+ T cells
    tcells <- subset(seu_obj, celltype == "T cells")
    tcells$cd3 <- colMeans(GetAssayData(tcells, slot = "data")[c("CD3D", "CD3E", "CD3G"), ]) > 0
    tcells$cd4 <- GetAssayData(tcells, slot = "data")[c("CD4"), ] > 0
    tcells$cd8 <- colMeans(GetAssayData(tcells, slot = "data")[c("CD8A", "CD8B"), ]) > 0
    table(tcells$cd4, tcells$cd8, tcells$cd3)
    tcells$sub_celltype <- ifelse((tcells$cd3 & tcells$cd8 & !(tcells$cd4)), "CD8_T",
                                   ifelse((tcells$cd3 & !(tcells$cd8) & tcells$cd4), "CD4_T", "Unknown_T"))
    print(table(tcells$sub_celltype))
    
    tmp <- tcells$sub_celltype[match(colnames(seu_obj), colnames(tcells))]
    seu_obj$celltype2 <- ifelse(seu_obj$celltype == "T cells", tmp, seu_obj$celltype)
    table(seu_obj$sub_celltype, seu_obj$celltype2, useNA = "ifany")
    
    ## count cells for cell type x sample
    res = as.data.frame.matrix(table(seu_obj$samID, seu_obj$celltype2))
    write.table(res, file = paste0("res/", dir.used, tumor, ".tsv"), quote = F, sep = "\t")
    
    saveRDS.pigz(seu_obj, file = paste0("data/rds/TpN/", tumor, ".TpN.rds"), n.cores = 10)
}