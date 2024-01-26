rm(list = ls())
library(Seurat)
library(scMetab)
library(RColorBrewer)
library(fastSave)
library(ggplot2)

seed.use = 2021
tmp.text <- element_text(family="sans", size=18)

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

seu.list <- lapply(1:10, function(i) {
    tumor <- dataset$DataSets[i]
    print(tumor)
    seu_obj <- readRDS(paste0("../", dataset$Directory[i], "/data/rds/", dataset$Prefix[i], ".TpN.merged.rds"))
    
    set.seed(seed.use)
    seu_obj <- seu_obj[, sample(x = colnames(seu_obj), size = 0.1*dim(seu_obj)[2])]
    seu_obj$dataset <- tumor
    
    gc()
    return(seu_obj)
})
set.seed(NULL)
names(seu.list) <- dataset$DataSets
nCells <- sapply(seu.list, function(x) dim(x)[2])
sum(nCells)

## clustering using highly variable genes
seu <- merge(seu.list[[1]], seu.list[-1])
seu <- seurat_pipe(seu, useHarmony = T, harmony.group.var = "dataset", seed.use = seed.use, 
                   fast_tsne_path = "D:/software/FIt-SNE-master/bin/FItSNE.exe")
saveRDS(seu, file = "data/rds/all_dataset_harmony.rds")

pdf(paste0("plot/1_clustering/All_dataset_merged/Dimplot_harmony.pdf"), width = 12, height = 10)
print(DimPlot(seu, reduction = "umap", group.by = "celltype", cols = celltype_color, 
              pt.size = 0.8, raster = T) +
          theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
print(DimPlot(seu, reduction = "tsne", group.by = "celltype", cols = celltype_color, 
              pt.size = 0.8, raster = T) +
          theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
print(DimPlot(seu, reduction = "umap", group.by = "dataset", cols = dataset_color, 
              pt.size = 0.8, raster = T) +
          theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
print(DimPlot(seu, reduction = "tsne", group.by = "dataset", cols = dataset_color, 
              pt.size = 0.8, raster = T) +
          theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
dev.off()

## clustering using metabolic genes
seu.metab <- merge(seu.list[[1]], seu.list[-1])
metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu.metab))
seu.metab <- seurat_pipe(seu.metab, features = metab.gene, 
                         useHarmony = T, harmony.group.var = "dataset", seed.use = seed.use, 
                         fast_tsne_path = "D:/software/FIt-SNE-master/bin/FItSNE.exe")
saveRDS(seu.metab, file = "data/rds/all_dataset_harmony_metab_features.rds")


pdf(paste0("plot/1_clustering/All_dataset_merged/Dimplot_harmony_metab_features.pdf"), 
    width = 12, height = 10)
print(DimPlot(seu.metab, reduction = "umap", group.by = "celltype", cols = celltype_color, 
              pt.size = 0.8, raster = T) +
          theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
print(DimPlot(seu.metab, reduction = "tsne", group.by = "celltype", cols = celltype_color, 
              pt.size = 0.8, raster = T) +
          theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
print(DimPlot(seu.metab, reduction = "umap", group.by = "dataset", cols = dataset_color, 
              pt.size = 0.8, raster = T) +
          theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
print(DimPlot(seu.metab, reduction = "tsne", group.by = "dataset", cols = dataset_color, 
              pt.size = 0.8, raster = T) +
          theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text))
dev.off()

## Plot top cell type-specific metabolic genes
top_genes <- read.csv("res/1_clustering/PAS_markers/Celltype_top_genes.csv")
DimPlot(seu.metab, group.by = "celltype", raster = T, reduction = "tsne", cols = celltype_color)
FeaturePlot(seu.metab, features = top_genes$x, raster = T, reduction = "tsne", 
            cols = rev(brewer.pal(11, "Spectral")))
ggsave(filename = "top_metabolic_genes.pdf", height = 24, width = 16)
