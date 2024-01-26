library(reshape2)

gene_markers.list <- lapply(rownames(dataset), function(tumor) {
    gene_markers <- read.delim(paste0("res/1_clustering/PAS_clustering/", tumor, "_Celltype_Gene_Markers_All.tsv"))
    gene_markers$log10padj <- (-1) * sign(gene_markers$avg_log2FC) * log10(gene_markers$p_val_adj)
    gene_markers$tumor <- tumor
    
    return(gene_markers)
})

names(gene_markers.list) <- rownames(dataset)
celltype_gene_markers <- do.call(rbind, gene_markers.list)

a <- subset(celltype_gene_markers, cluster == "Malignant")
b <- reshape2::dcast(a, gene ~ tumor, value.var = "avg_log2FC")




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
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

seu.list <- lapply(1:nrow(dataset), function(i) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(paste0("data/rds/PAS/", tumor, ".TpN.PAS.unintegrated.rds"))
    seu_obj <- subset(seu_obj, TorN == "T")
    return(seu_obj)
})
names(seu.list) <- rownames(dataset)


seu.list <- seu.list[sort(names(seu.list))]

## plot function
plotfun <- function(feature, assay, plot_box = FALSE) {
    tmp.text <- element_text(family="sans", size=12)
    p = lapply(names(seu.list), function(x) {
        seu <- seu.list[[x]]
        DefaultAssay(seu) <- assay
        tmpp = VlnPlot(seu, features = feature, group.by = "celltype", 
                pt.size = 0, cols = celltype_color, combine = F)[[1]] +
            ggtitle(x) + NoLegend() + 
            ylab("") + xlab("") +
            theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, axis.title = tmp.text)
        if (plot_box) {
            tmpp = tmpp + 
                geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9))
        }
        return(tmpp)
    })
    return(wrap_plots(p, ncol = 2) + plot_annotation(feature))
}

## Plot VlnPlot for SOD3 and LDHB
p1.list <- plotfun(feature = "SOD3", assay = "RNA")
p2.list <- plotfun(feature = "LDHB", assay = "RNA")

pdf("plot/1_clustering/SOD3_LDHB_vln.pdf", width = 12, height = 12)
wrap_plots(p1, p2, ncol = 2)
dev.off()

## plot VlnPlot for GLS, GLS2, GLUD1, GLUL
p1.list <- lapply(c("Glycolysis and Gluconeogenesis", "Oxidative Phosphorylation", 
                    "Glutamate metabolism"), plotfun, assay = "Metab", plot_box = TRUE)
p2.list <- lapply(c("GLS", "GLS2", "GLUD1", "GLUL"), plotfun, assay = "RNA")
pdf("plot/1_clustering/glutamate_vln.pdf", width = 6, height = 12)
p1.list
p2.list
dev.off()


# p2.list <- lapply(seu.list, function(seu) {
#     
#     DimPlot(seu, group.by = "celltype", reduction = "tsne", raster = T, 
#             pt.size = 0.8, cols = celltype_color) + NoLegend() 
# })
# wrap_plots(p2.list, ncol = 2)
# ggsave("celltype_tsne.pdf", width = 5, height = 12)
