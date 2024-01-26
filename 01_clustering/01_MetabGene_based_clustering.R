# ==================================================================================================
# Author: Zhou Zhe
# Program: Metabolic gene based clustering
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
seed.use = 2021
dir.used <- "1_clustering/MetabGene_clustering/"
dir.create(file.path("res", dir.used), recursive = T)
dir.create(file.path("plot", dir.used), recursive = T)
# file.prefix <- ""
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
used_celltype = c("Malignant", "Epithelial", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")

# skip i = 1 !!!
for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(paste0("data/rds/TpN/", tumor, ".TpN.rds"))
    metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu_obj))
    if (!("percent.metab" %in% colnames(seu_obj[[]]))) {
        print("No percent.metab!!")
        seu_obj$percent.metab <- PercentageFeatureSet(seu_obj, features = metab.gene)
    }
    
    
    # plot nFeature_metab
    tmpfile <- paste0("plot/", dir.used, tumor, "_nFeature_metab.pdf")
    tmp_seu <- subset(seu_obj, celltype %in% used_celltype)
    tmp.text <- element_text(family="sans", size=18)
    pdf(tmpfile, width = 12, height = 4)
    print(VlnPlot(tmp_seu, features = c("nFeature_metab", "nCount_metab", "percent.metab"),
            group.by = "celltype", cols = celltype_color, pt.size = 0) &
        geom_boxplot(width = 0.1, outlier.shape = NA) &
        theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(VlnPlot(tmp_seu, features = c("nFeature_metab", "nCount_metab", "percent.metab"),
            group.by = "celltype", split.by = "TorN", cols = TvN_color, pt.size = 0) &
        geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
        theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    dev.off()
    
    # II. Clustering for Metab Gene ========
    seu_obj <- seurat_pipe(seu_obj, features = metab.gene, seed.use = seed.use)
    VariableFeatures(seu_obj) <- metab.gene
    saveRDS.pigz(seu_obj, file = paste0("data/rds/MetabGene/", tumor, ".TpN.unintegrated.rds"), n.cores = 10)

    tmp.text <- element_text(family="sans", size=18)
    pdf(paste0("plot/", dir.used, tumor, "_MetabGene_unintegrated_Imp.pdf"), width = 12, height = 10)
    print(DimPlot(seu_obj, reduction = "umap", group.by = "celltype", cols = celltype_color, pt.size = 0.2) +
              theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "umap", group.by = "celltype", cols = celltype_color, pt.size = 0.8, raster = T) +
              theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "tsne", group.by = "celltype", cols = celltype_color, pt.size = 0.2) +
              theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "tsne", group.by = "celltype", cols = celltype_color, pt.size = 0.8, raster = T) +
              theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    dev.off()
}
