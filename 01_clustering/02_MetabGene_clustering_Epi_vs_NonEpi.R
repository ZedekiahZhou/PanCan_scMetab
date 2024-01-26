# ==================================================================================================
# Author: Zhou Zhe
# Program: Metabolic gene based clustering, Epi and nonEpi splited for BRCA
# Version: 1.0
# Date: Mar 5, 2023
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
seed.use = 2021
dir.used <- "1_clustering/MetabGene_clustering/"
dir.create(file.path("res", dir.used), recursive = T)
dir.create(file.path("plot", dir.used), recursive = T)
# file.prefix <- ""
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
used_celltype = c("Malignant", "Epithelial", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")

# skip i = 1 !!!
for (i in 1:1) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(paste0("../", dataset$Directory[i], "/data/rds/", dataset$Prefix[i], ".TpN.merged.rds"))
    metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu_obj))
    # seu_obj <- subset(seu_obj, TorN == "T")
    
    malig <- subset(seu_obj, celltype == "Malignant" & TorN == "T")
    nonmalig <- subset(seu_obj, celltype != "Malignant" & TorN == "T")
    
    # II. Clustering for Metab Gene ========
    malig <- seurat_pipe(malig, features = metab.gene, seed.use = seed.use, 
                         fast_tsne_path = "D:/software/FIt-SNE-master/bin/FItSNE.exe")
    nonmalig <- seurat_pipe(nonmalig, features = metab.gene, seed.use = seed.use, 
                         fast_tsne_path = "D:/software/FIt-SNE-master/bin/FItSNE.exe")
    
    tumorType_color <- celltype_color[c(1, 5, 6)]
    names(tumorType_color) <- c("TNBC", "ER_BC", "HER2_BC")
    patient_color <- colorRampPalette(celltype_color)(44)
    names(patient_color) <- sort(unique(seu_obj$patientID))
    saveRDS(patient_color, file = "data/rds/patient_color/BRCA_patient_color.rds")
    patient_color <- patient_color[names(patient_color) %in% malig$patientID]
    
    tmp.text <- element_text(family="sans", size=6)
    pdf(paste0("plot/", dir.used, "BRCA_EPi_vs_NonEpi.pdf"), width = 9.6, height = 3.2)
    p1 = DimPlot(malig, cols = celltype_color, reduction = "tsne", group.by = "celltype", 
                 pt.size = 0.8, raster = T) + 
        DimPlot(malig, cols = tumorType_color, reduction = "tsne", group.by = "tumorType", 
                pt.size = 0.8, raster = T) +
        DimPlot(malig, cols = patient_color, reduction = "tsne", group.by = "patientID", 
                pt.size = 0.8, raster = T) & 
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text, 
              legend.key.size = unit(1, "mm"))
    print(p1)
    p2 = DimPlot(nonmalig, cols = celltype_color, reduction = "tsne", group.by = "celltype", 
                 pt.size = 0.8, raster = T) + 
        DimPlot(nonmalig, cols = tumorType_color, reduction = "tsne", group.by = "tumorType", 
                pt.size = 0.8, raster = T) +
        DimPlot(nonmalig, cols = patient_color, reduction = "tsne", group.by = "patientID", 
                pt.size = 0.8, raster = T) & 
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text, 
              legend.key.size = unit(1, "mm"))
    print(p2)
    dev.off()
        
    saveRDS(malig, file = paste0("data/rds/MetabGene/BRCA_malig_unintegrated.rds"))
    saveRDS(nonmalig, file = paste0("data/rds/MetabGene/BRCA_nonmalig_unintegrated.rds"))
    
}