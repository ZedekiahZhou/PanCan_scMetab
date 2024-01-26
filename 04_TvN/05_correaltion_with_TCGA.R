# ==================================================================================================
# Author: Zhou Zhe
# Program: compare celltype TvN to TCGA TvN
#          Datasets: Cell_2021_Pelka_colon
# Version: 1.0
# Date: Feb 14, 2022
# ==================================================================================================

rm(list=ls())
library(Seurat)
library(scMetab)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ComplexHeatmap)

seed.use = 2021
dir.used <- "2_TvN/TCGA/"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
seed.use = 2021

# I. Get All log2FC matrix -----------------------------------------------------------------------------------
log2FC_list <- lapply(c(1, 4, 7:10), function(i) {
    tumor = dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    seu_obj <- readRDS(paste0("../", dataset$Directory[i], "/data/rds/", dataset$Prefix[i], ".TpN.merged.rds"))
    if (!("samID" %in% colnames(seu_obj[[]]))) seu_obj$samID = seu_obj$patientID
    if (tumor == "PRAD") seu_obj$samID <- seu_obj$samID2
    
    TCGA_res <- read.delim(paste0("res/2_TvN/TCGA/", tumor, "_Gene_markers.tsv"))
    rownames(TCGA_res) <- TCGA_res$symbol
    
    metab.gene <- intersect(rownames(TCGA_res), gaudeMetabDf$GeneSymbol)
    metab.gene <- intersect(rownames(seu_obj), metab.gene)
    
    
    ###### ====== correlation to celltype TvN ====== ######
    celltypes <- c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")
    gene_stat <- sapply(celltypes, function(x) {
        tmpfile <- paste0("res/2_TvN/AUCell/", tumor, "/", tumor, "_", sub(" ", "_", x), "_Gene_Markers.tsv")
        gene_markers <- read.delim(file = tmpfile)
        rownames(gene_markers) <- gene_markers$gene
        gene_markers <- gene_markers[metab.gene, ]
        return(gene_markers[, "avg_log2FC"])
    })
    rownames(gene_stat) <- metab.gene
    gene_stat <- as.data.frame(gene_stat)
    
    
    ######## ====== pseudo bulk ====== ########
    seu_obj <- NormalizeData(seu_obj)
    pseudo_sam <- AggregateExpression(seu_obj, slot = "data", features = metab.gene, 
                                      group.by = "samID", return.seurat = T)
    pseudo_sam$TorN <- seu_obj$TorN[match(colnames(pseudo_sam), seu_obj$samID)]
    Idents(pseudo_sam) <- "TorN"
    pseudo_markers <- FindMarkers(pseudo_sam, slot = "data", assay = "RNA",
                                  ident.1 = "T", ident.2 = "N", 
                                  logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                                  max.cells.per.ident = Inf, random.seed = seed.use)
    
    
    ######## ====== Pooled all log2FC ====== ########
    stat_matrix <- cbind(gene_stat[metab.gene, ], pseudo_bulk = pseudo_markers[metab.gene, "avg_log2FC"], 
                         TCGA = TCGA_res[metab.gene, "log2FC"])
    
    stat_matrix <- stat_matrix[!is.na(stat_matrix$TCGA), ]
    # cor(stat_matrix, method = "spearman")
    # cor(stat_matrix, method = "pearson")
    print(colSums(stat_matrix != 0))
    
    write.csv(stat_matrix, file = paste0("res/", dir.used, tumor, "_celltype_and_TCGA_log2FC.csv"))
    return(stat_matrix)
})
names(log2FC_list) <- dataset$DataSets[c(1, 4, 7:10)]



####### ================ TCGA Corrrelation ==================== ###################
celltypes <- c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts", "pseudo_bulk")
mat <- sapply(log2FC_list, function(x) {
    cor(x, method = "spearman")[celltypes, "TCGA"]
})

tmpfile <- paste0("plot/", dir.used, "PanCan_correlation_with_TCGA.pdf")
tmp.text <- gpar(fontfamily="sans", fontsize = 8)
myCols <- brewer.pal(n=11, name="RdBu")[c(8, 6, 2)]
pdf(tmpfile, width = 5, height = 1.5)
Heatmap(mat, circlize::colorRamp2(c(-0.1, 0, 0.5), myCols), 
        rect_gp = gpar(col = "white"), 
        row_names_max_width = unit(6, "cm"),
        row_names_gp = tmp.text, column_names_gp = tmp.text, 
        column_names_rot = 45, 
        row_title_gp = tmp.text, column_title_gp = tmp.text, 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        name = "Correlatino\nwith TCGA", 
        cluster_rows = T, cluster_columns = FALSE)
dev.off()

ggscatter(log2FC_list[[1]], x = "pseudo_bulk", y = "TCGA",
          add = "reg.line", conf.int = TRUE,
          add.params = list(fill = "red3")) +
    stat_cor(method = "spearman")












