# ==================================================================================================
# Author: Zhou Zhe
# Program: prepare gmt file of celltype signatures for SCENIC AUCell
#          Datasets: for all datasets
#          Samples: ..
#          TorN: only Tumor samples
#          Celltype: All celltype included
# Version: 1.0
# Date: Aug 2, 2022
# ==================================================================================================

# init -------
rm(list=ls())


library(tidyverse)
library(Seurat)
library(SCopeLoomR)
library(scMetab)
library(GSEABase)

dir.used <- "res/5_SCENIC/"
file.prefix <- "SCENIC_prepare_"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

# signature to gmt
sig2gmt <- function(sig.df, bg.genes, outfile) {
    sig.list = split(sig.df$GeneSymbol, f = sig.df$Pathway)
    filtered.sig.list = lapply(sig.list, function(x) intersect(x, bg.genes))
    
    rawLen = sapply(sig.list, length)
    filteredLen = sapply(filtered.sig.list, length)
    missPath = names(rawLen)[filteredLen/rawLen <= 0.8]
    cat("**NOTE** Following pathways have less than 80% signatures present in data:\n", missPath)
    
    sig.gmt <- lapply(names(filtered.sig.list), function(x) {
        GeneSet(filtered.sig.list[[x]], setName = x)
    })
    sig.gmt <- GeneSetCollection(sig.gmt)
    toGmt(sig.gmt, outfile)
    return(sig.gmt)
}



for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    sink(paste0(dir.used, "log/", file.prefix, tumor, ".log"))
    
    loom = open_loom(file.path = paste0(dir.used, "whole_loom/", tumor, "_filtered_whole.loom"))
    genes = get_genes(loom)
    close_loom(loom)
        
    # filter gene set gmt file
    a = sig2gmt(sig.df = tcellSig, bg.genes = genes, outfile = paste0(dir.used, "dataset_gmt/", tumor, "_tcellSig.gmt"))
    b = sig2gmt(sig.df = myeloidSig, bg.genes = genes, outfile = paste0(dir.used, "dataset_gmt/", tumor, "_myeloidSig.gmt"))
    c = sig2gmt(sig.df = cancerHallmarksDf, bg.genes = genes, outfile = paste0(dir.used, "dataset_gmt/", tumor, "_cancerHallmarks.gmt"))
    toGmt(GeneSetCollection(c(a, b, c)), paste0(dir.used, "dataset_gmt/", tumor, "_celltype_sig.gmt"))
    sink()
}

