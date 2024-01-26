# ==================================================================================================
# Author: Zhou Zhe
# Program: prepare whole data (include all genes, filter low expressed) and filtered data (only TF
#          and metabolic genes) for SCENIC
#          Datasets: for all datasets
#          Samples: ..
#          TorN: Tumor and Normal samples
#          Celltype: All celltype included
# Version: 1.0
# Date: Jul 7, 2022
# ==================================================================================================

# init -------
rm(list=ls())


library(tidyverse)
library(Seurat)
library(SCopeLoomR)
library(scMetab)
library(GSEABase)

dir.used <- "res/02_SCENIC/"
file.prefix <- "SCENIC_prepare_"

########### !!!!!!!!
for (folder in c("log", "metab_loom", "whole_loom", "malignant_loom", "T_loom", "myeloid_loom")) {
    dir.create(file.path(dir.used, folder), recursive = T)
}

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

# metabolic genes
metab.genes <- unique(gaudeMetabDf$GeneSymbol)
metab.list <- split(gaudeMetabDf$GeneSymbol, f = gaudeMetabDf$Pathway)
metab.len <- sapply(metab.list, length)
metab.list <- metab.list[metab.len >=5]

# write out metabolic genes to GMT file
metab.gmt <- lapply(names(metab.list), function(x) {
    GeneSet(metab.list[[x]], setName = x)
})
metab.gmt <- GeneSetCollection(metab.gmt)
toGmt(metab.gmt, "data/GeneSets/gaudeMetab.gmt")


for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    sink(paste0(dir.used, "log/", file.prefix, tumor, ".log"))
    
    seu_obj <- readRDS(paste0("data/rds/TpN/", tumor, ".TpN.rds"))
    # seu_obj <- readRDS(paste0("../", dataset$Directory[i], "/data/rds/", dataset$Prefix[i], ".TpN.merged.rds"))
    cat(paste0("There are ", dim(seu_obj)[1], " genes and ", dim(seu_obj)[2], " cells in ", 
               dataset$Prefix[i], ".TpN.merged.rds\n\n"))
    
    
    # filter low expressed genes
    cat("Filter low expressed genes ====== \n")
    percents <- ClusterPercent(seu_obj, assay = "RNA", group.by = "celltype")
    max_percent <- sort(apply(percents, 1, max), decreasing = TRUE)
    cat("Quantile of gene's max expressed percent across celltypes: \n")
    print(quantile(max_percent, seq(0, 1, 0.05)))
    cat(sum(max_percent > 0.01), "genes expressed in more than 1% cells of at least one celltype!\n\n")
    
    # save whole loom
    max_percent <- max_percent[max_percent > 0.01]
    seu_obj <- seu_obj[names(max_percent), ]
    build_loom(file.name = paste0(dir.used, "whole_loom/", tumor, "_filtered_whole.loom"),
               dgem = counts,
               title = dataset$DataSets[i])
    
    
    # filter metabolic genes
    filtered.metab.list = lapply(metab.list, function(x) intersect(x, names(max_percent)))
    metab.gmt <- lapply(names(filtered.metab.list), function(x) {
        GeneSet(filtered.metab.list[[x]], setName = x)
    })
    metab.gmt <- GeneSetCollection(metab.gmt)
    toGmt(metab.gmt, paste0(dir.used, "dataset_gmt/", tumor, "_gaudeMetab.gmt"))
    
    
    # filter TF and metab genes 
    cat("Filter TF and metab genes ======== \n")
    cat(paste0("There are ", length(TFgenes), " TF genes, and ", 
               length(intersect(rownames(seu_obj), TFgenes)),  " of them are in the dataset.\n"))
    cat(paste0("There are ", length(unique(gaudeMetabDf$GeneSymbol)), " metab genes, and ", 
               length(intersect(rownames(seu_obj), unique(gaudeMetabDf$GeneSymbol))),  " of them are in the dataset.\n"))
    seu_obj <- seu_obj[rownames(seu_obj) %in% c(gaudeMetabDf$GeneSymbol, TFgenes), ]
    cat(paste0("There are ", dim(seu_obj)[1], " genes in the dataset which are either TF or metab genes.\n\n"))
    
    
    # expressed features of each cell
    counts <- GetAssayData(seu_obj, slot = "counts")
    nFeature <- colSums(counts > 0)
    cat("Quantile of cell's nFeature: \n")
    print(quantile(nFeature, seq(0, 1, 0.05)))
    
    # save TF and metab genes loom
    build_loom(file.name = paste0(dir.used, "metab_loom/", tumor, "_filtered_scenic.loom"), 
               dgem = counts, 
               title = dataset$DataSets[i])
    
    sink()
}

# for (tumor in rownames(dataset)) {
#     system(paste("scp", file.path("/store/Pancancer_Metab/Public", dataset[tumor, ]$dir, dir.used, 
#                                   paste0(tumor, "_filtered_scenic.loom")), 
#                  file.path("shengke1:~/Proj/pySCENIC/data", dataset[tumor, ]$dir)))
# }
