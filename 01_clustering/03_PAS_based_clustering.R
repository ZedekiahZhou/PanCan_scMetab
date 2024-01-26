# ==================================================================================================
# Author: Zhou Zhe
# Program: PAS clustering and celltype specific PAS
# Version: 1.0
# Date: Jul 13, 2022
# ==================================================================================================

# 0. INIT ========
rm(list=ls())

library(dplyr)
library(tidyverse)
library(tidyselect)
library(Seurat)
library(scMetab)
library(RColorBrewer)
library(SCENIC)
library(data.table)
library(fastSave)
library(patchwork)

n_cores <- 4
seed.use = 2021
dir.used <- "1_clustering/PAS_clustering/"
# file.prefix <- ""
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets


# skip i = 1 !!!
for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    # I. Load Data ========
    metab_pas <- data.frame(fread(paste0("res/5_SCENIC/aucell/", tumor, "_metab_aucell.csv")),
                            row.names = 1, check.names = F)
    TF_pas <- data.frame(fread(paste0("res/5_SCENIC/aucell/", tumor, "_TF_aucell.csv")),
                         row.names = 1, check.names = F)
    seu_obj <- readRDS(paste0("data/rds/TpN/", tumor, ".TpN.rds"))
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj <- NormalizeData(seu_obj)
    seu_obj$celltype_x_TorN <- paste(seu_obj$celltype, seu_obj$TorN, sep = "_")

    seu_obj[["Metab"]] <- CreateAssayObject(counts = t(metab_pas))
    seu_obj[["TF"]] <- CreateAssayObject(counts = t(TF_pas))

    # II. Clustering for metab ========
    # cells = sample(colnames(seu_obj), size = 2000)
    DefaultAssay(seu_obj) <- "Metab"
    seu_obj <- ScaleData(seu_obj)
    VariableFeatures(seu_obj) <- rownames(seu_obj)
    seu_obj <- seu_obj %>%
        RunPCA(seed.use = seed.use, verbose = F) %>%
        RunUMAP(dims = 1:30, seed.use = seed.use, verbose = F) %>%
        RunTSNE(dims = 1:30, seed.use = seed.use, tsne.method = "FIt-SNE", nthreads = 10, verbose = F)
    gc()

    tmp.text <- element_text(family="sans", size=18)
    pdf(paste0("plot/", dir.used, tumor, "_Metab_PAS_unintegrated_Imp.pdf"), width = 12, height = 10)
    print(DimPlot(seu_obj, reduction = "umap", group.by = "celltype", cols = celltype_color, pt.size = 0.2) +
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "umap", group.by = "celltype", cols = celltype_color, pt.size = 0.8, raster = T) +
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "tsne", group.by = "celltype", cols = celltype_color, pt.size = 0.2) +
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "tsne", group.by = "celltype", cols = celltype_color, pt.size = 0.8, raster = T) +
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    dev.off()


    # III. Clustering for TF ========
    # cells = sample(colnames(seu_obj), size = 2000)
    DefaultAssay(seu_obj) <- "TF"
    seu_obj <- ScaleData(seu_obj)
    VariableFeatures(seu_obj) <- rownames(seu_obj)
    seu_obj <- seu_obj %>%
        RunPCA(seed.use = seed.use, verbose = F) %>%
        RunUMAP(dims = 1:30, seed.use = seed.use, verbose = F) %>%
        RunTSNE(dims = 1:30, seed.use = seed.use, tsne.method = "FIt-SNE", nthreads = 10, verbose = F)
    gc()

    tmp.text <- element_text(family="sans", size=18)
    pdf(paste0("plot/", dir.used, tumor, "_TF_PAS_unintegrated_Imp.pdf"), width = 12, height = 10)
    print(DimPlot(seu_obj, reduction = "umap", group.by = "celltype", cols = celltype_color, pt.size = 0.2) +
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "umap", group.by = "celltype", cols = celltype_color, pt.size = 0.8, raster = T) +
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "tsne", group.by = "celltype", cols = celltype_color, pt.size = 0.2) +
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    print(DimPlot(seu_obj, reduction = "tsne", group.by = "celltype", cols = celltype_color, pt.size = 0.8, raster = T) +
        theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
    dev.off()

    saveRDS.pigz(seu_obj, file = paste0("data/rds/PAS/", tumor, ".TpN.PAS.unintegrated.rds"), n.cores = 10)


    # IV. celltype markers ======
    seu_tumor <- subset(seu_obj, TorN == "T")
    Idents(seu_tumor) <- "celltype"
    metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu_tumor@assays$RNA))
    metab_markers <- FindAllMarkers(seu_tumor, slot = "data", assay = "Metab",
                                  logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                                  return.thresh = 2, 
                                  max.cells.per.ident = 1000, random.seed = seed.use)
    TF_markers <- FindAllMarkers(seu_tumor, slot = "data", assay = "TF",
                                 logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                                 return.thresh = 2, 
                                 max.cells.per.ident = 1000, random.seed = seed.use)
    gene_markers <- FindAllMarkers(seu_tumor, slot = "data", assay = "RNA", features = metab.gene,
                                   logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                                   return.thresh = 2, 
                                   max.cells.per.ident = 1000, random.seed = seed.use)
    tmpfile <- paste0("res/", dir.used, tumor)
    write.table(metab_markers, file = paste0(tmpfile, "_Celltype_Metab_Markers_All.tsv"), quote = F, row.names = F, sep = "\t")
    write.table(TF_markers, file = paste0(tmpfile, "_Celltype_TF_Markers_All.tsv"), quote = F, row.names = F, sep = "\t")
    write.table(gene_markers, file = paste0(tmpfile, "_Celltype_Gene_Markers_All.tsv"), quote = F, row.names = F, sep = "\t")

    # V. RSS ========
    pdf(paste0("plot/", dir.used, tumor, "_RSS.pdf"), width = 12, height = 10)
    metab_pas <- metab_pas[colnames(seu_tumor), ]
    metab_rss <- calcRSS(t(metab_pas), cellAnnotation = as.factor(seu_tumor$celltype))
    print(plotRSS(metab_rss))
    p = lapply(unique(seu_tumor$celltype), function(x) plotRSS_oneSet(metab_rss, setName = x))
    print(Reduce("+", p))
    
    TF_pas <- TF_pas[colnames(seu_tumor), ]
    TF_rss <- calcRSS(t(TF_pas), cellAnnotation = as.factor(seu_tumor$celltype))
    print(plotRSS(TF_rss))
    p = lapply(unique(seu_tumor$celltype), function(x) plotRSS_oneSet(TF_rss, setName = x))
    print(Reduce("+", p))
    dev.off()
}    