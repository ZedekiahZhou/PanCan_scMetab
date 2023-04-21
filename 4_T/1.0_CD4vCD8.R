# ====================================================================
# Author: Zhou Zhe
# Function: pathway correlation of CD8_T
# Version: 1.0
# Date: Aug 3, 2022
# ====================================================================

rm(list=ls())
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)
library(data.table)
library(scMetab)
library(GSEABase)
library(tidyr)

pseudo <- 1e-50
seed.use = 2021

dir.used <- "4_T/CD4vCD8/"
dir.create(file.path("res", dir.used), recursive = T)
dir.create(file.path("plot", dir.used), recursive = T)

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

# I. CD4vCD8 DE -----------------------------------------------------------------------------------
metab_marker_list <- lapply(1:nrow(dataset), function(i) {
    tumor = dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    seu_obj <- readRDS(paste0("data/rds/PAS/", tumor, ".TpN.PAS.unintegrated.rds"))
    print(table(seu_obj$celltype, seu_obj$TorN))
    seu_obj <- subset(seu_obj, celltype == "T cells" & TorN == "T")
    
    ## CD3, CD4, CD8
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj$cd3 <- colMeans(GetAssayData(seu_obj, slot = "data")[c("CD3D", "CD3E", "CD3G"), ]) > 0
    seu_obj$cd4 <- GetAssayData(seu_obj, slot = "data")[c("CD4"), ] > 0
    seu_obj$cd8 <- colMeans(GetAssayData(seu_obj, slot = "data")[c("CD8A", "CD8B"), ]) > 0
    table(seu_obj$cd4, seu_obj$cd8, seu_obj$cd3)
    seu_obj$sub_celltype <- ifelse((seu_obj$cd3 & seu_obj$cd8 & !(seu_obj$cd4)), "CD8_T",
                                   ifelse((seu_obj$cd3 & !(seu_obj$cd8) & seu_obj$cd4), "CD4_T", "Unknown_T"))
    print(table(seu_obj$sub_celltype))
    seu_obj <- subset(seu_obj, sub_celltype != "Unknown_T")
    
    # DE
    DefaultAssay(seu_obj) <- "Metab"
    Idents(seu_obj) <- "sub_celltype"
    
    # metab markers
    metab_markers <- FindMarkers(seu_obj, slot = "data", assay = "Metab",
                                 ident.1 = "CD4_T", ident.2 = "CD8_T",
                                 logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                                 random.seed = seed.use)
    metab_markers$pathway <- rownames(metab_markers)
    
    # TF markers
    TF_markers <- FindMarkers(seu_obj, slot = "data", assay = "TF",
                              ident.1 = "CD4_T", ident.2 = "CD8_T",
                              logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                              random.seed = seed.use)
    TF_markers$pathway <- rownames(TF_markers)
    
    
    # Metabolic gene markers
    # gene_markers <- FindMarkers(seu_obj, slot = "data", assay = "RNA", features = rownames(seu_obj[["RNA"]]),
    #                             ident.1 = "CD4_T", ident.2 = "CD8_T",
    #                             logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
    #                             max.cells.per.ident = 1000, random.seed = seed.use)
    # gene_markers$gene <- rownames(gene_markers)
    
    # save data
    write.table(metab_markers, file = paste0("res/", dir.used, tumor, "CD4vCD8_Metab_Markers.tsv"), quote = F, row.names = F, sep = "\t")
    write.table(TF_markers, file = paste0("res/", dir.used, tumor, "CD4vCD8_TF_Markers.tsv"), quote = F, row.names = F, sep = "\t")
    # write.table(gene_markers, file = paste0("res/", dir.used, tumor, "_CD4vCD8_Gene_Markers.tsv"), quote = F, row.names = F, sep = "\t")
    return(metab_markers)
})
names(metab_marker_list) <- dataset$DataSets


# 2. Metabolic Markers Plot -----------------------------------------------------------------------
metab_markers <- lapply(c(1:10), function(i) {
    tumor = dataset$DataSets[i]
    tmp.data <- read.delim(paste0("res/4_T/CD4vCD8/", tumor, "CD4vCD8_Metab_Markers.tsv"))
    tmp.data$tumor <- tumor
    tmp.data
})
metab_markers <- do.call(rbind, metab_markers)
metab_markers <- subset(metab_markers, (p_val_adj < 0.01) & (abs(avg_log2FC) > 0.01) & ((pct.1 > 0.1) | (pct.2 > 0.1)))

metab_markers$log10padj <- sign(metab_markers$avg_log2FC) * (-1) * log10(metab_markers$p_val_adj)
metab_markers$log10padj_pseudo <- sign(metab_markers$avg_log2FC) * (-1) * log10(metab_markers$p_val_adj + pseudo)
metab_markers <- metab_markers[order(metab_markers$avg_log2FC, decreasing = T), ]
metab_markers$pathway <- factor(metab_markers$pathway, levels = unique(metab_markers$pathway))

tmpfile <- paste0("plot/", dir.used, "CD4vCD8_Metab_dotplot.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 4.8, height = ceiling(length(unique(metab_markers$pathway))/8))
print(ggplot(metab_markers, aes(x = tumor, y = pathway, color = log10padj_pseudo)) + 
          geom_point(aes(size = abs(avg_log2FC))) + 
          scale_color_distiller(palette = "RdBu") + 
          scale_size(limits = c(0, 0.5), range = c(0.5, 3.5)) + 
          theme_classic() + 
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black")))
dev.off()
