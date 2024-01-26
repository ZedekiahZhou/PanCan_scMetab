# ==================================================================================================
# Author: Zhou Zhe
# Program: all gene based clustering
# Version: 1.0
# Date: Aug 12, 2022
# ==================================================================================================

# 0. INIT ========
rm(list=ls())

library(tidyverse)
library(Seurat)
library(scMetab)
library(RColorBrewer)
library(data.table)
library(fastSave)
library(patchwork)

n_cores <- 4
seed.use = 20201
dir.used <- "00_reproduce/AllGene_clustering/"
dir.create(file.path("res", dir.used), recursive = T)
dir.create(file.path("plot", dir.used), recursive = T)
# file.prefix <- ""
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets


for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(paste0("data/rds/TpN/", tumor, ".TpN.rds"))
    if (!("samID" %in% colnames(seu_obj[[]]))) seu_obj$samID = seu_obj$patientID
    if (tumor == "PRAD") seu_obj$samID <- seu_obj$samID2
    
    nSam = as.data.frame.array(table(seu_obj$samID, seu_obj$TorN))
    write.csv(nSam, file = paste0("res/", dir.used, tumor, "_nSam_nCell_of_TandN.csv"))
    
    # II. Clustering for All Gene ========
    seu_obj <- seurat_pipe(seu_obj, seed.use = seed.use)
    saveRDS.pigz(seu_obj, file = paste0("data/rds/AllGene/", tumor, ".TpN.unintegrated.rds"), n.cores = 10)
    
    
    
    tmp.text <- element_text(family="sans", size=18)
    pdf(paste0("plot/", dir.used, tumor, "_AllGene_unintegrated_Imp.pdf"), width = 12, height = 10)
    print(DimPlot(seu_obj, reduction = "umap", group.by = "celltype", cols = celltype_color, pt.size = 0.2) +
              theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "umap", group.by = "celltype", cols = celltype_color, pt.size = 0.8, raster = T) +
              theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "tsne", group.by = "celltype", cols = celltype_color, pt.size = 0.2) +
              theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "tsne", group.by = "celltype", cols = celltype_color, pt.size = 0.8, raster = T) +
              theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    dev.off()
    
    
    # Cell Type Propotion
    celltype_counts <- as.data.frame(table(seu_obj$celltype, seu_obj$TorN))
    colnames(celltype_counts) <- c("CellType", "Group", "Freq")
    celltype_counts$CellType <- factor(celltype_counts$CellType, levels = names(celltype_color))
    write.csv(celltype_counts, file = paste0("res/", dir.used, tumor, "_celltype_counts.csv"),
              row.names = F, quote = F)
}