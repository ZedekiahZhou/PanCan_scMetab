# ====================================================================
# Author: Zhou Zhe
# Function: pathway correlation of Myeloid
# Version: 1.0
# Date: Aug 3, 2022
# ====================================================================

rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)
library(data.table)
library(scMetab)
library(GSEABase)
library(tidyr)

pseudo <- 1e-50

dir.used <- "4_T/Mye/"
file.prefix <- "MyeloidSig_Cor_"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
myeloid.path <- unique(myeloidSig$Pathway)

pathway_cor_list <- lapply(1:nrow(dataset), function(i) {
    tumor = dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    seu_obj <- readRDS(paste0("../", dataset$Directory[i], "/data/rds/", dataset$Prefix[i], ".TpN.merged.rds"))
    table(seu_obj$celltype, seu_obj$TorN)    
    seu_obj <- subset(seu_obj, TorN == "T" & celltype == "Myeloid")
    
    # I. Load Data ========
    metab_pas <- data.frame(fread(paste0("res/5_SCENIC/aucell/", tumor, "_metab_aucell.csv")),
                            row.names = 1, check.names = F)[colnames(seu_obj), ]
    sig_pas <- data.frame(fread(paste0("res/5_SCENIC/aucell/", tumor, "_celltypeSig_aucell.csv")),
                          row.names = 1, check.names = F)[colnames(seu_obj), ]
    sig_pas <- sig_pas[, myeloid.path]
    gene.list <- c(geneIds(getGmt(paste0("res/5_SCENIC/dataset_gmt/", tumor, "_gaudeMetab.gmt"))), 
                   geneIds(getGmt(paste0("res/5_SCENIC/dataset_gmt/", tumor, "_myeloidSig.gmt"))))
    
    pathway_cor <- pathway_correlation(df = cbind(metab_pas, sig_pas), glist = gene.list)
    fastSave::saveRDS.pigz(pathway_cor, file = paste0("res/", dir.used, tumor, "Myeloid_pathway_correlation.rds"))
    return(pathway_cor)
})
names(pathway_cor_list) <- dataset$DataSets

pathway_cor_list <- lapply(dataset$DataSets, function(tumor) {
    return(readRDS(paste0("res/4_T/", tumor, "Myeloid_pathway_correlation.rds")))
})
names(pathway_cor_list) <- dataset$DataSets

get_pooled_data <- function(var) {
    res <- lapply(names(pathway_cor_list), function(tumor) {
        x <- pathway_cor_list[[tumor]][[var]]
        x <- data.frame(x[setdiff(rownames(x), myeloid.path), myeloid.path], check.names = FALSE)
        colnames(x) <- paste(colnames(x), tumor, sep = "_")
        return(x)
    })
    res <- Reduce(function(x, y) {
        z <- merge(x, y, by = "row.names", all = TRUE)
        rownames(z) <- z$Row.names
        z <- z[-1]
    } , res)
    return(res)
}
radj <- get_pooled_data("radj")
padj <- get_pooled_data("padj")

# keep correlation > 0.1 and adjust p < 0.05
idx <- (padj > 0.05) | (abs(radj) < 0.1) | (is.na(radj))
radj[idx] <- 0
radj <- radj[rowSums(radj != 0) > 0, ]
colnames(radj) <- sub("TIM_", "", colnames(radj))

# keep pathway have conserved correlation in >=5 tumor with at least 1 Signature
idx_df <- reshape2::melt(as.matrix(radj))
idx_df$Pathway <- sub("_.+", "", idx_df$Var2)
idx_df$Dataset <- sub("*[^_]+_", "", idx_df$Var2)
idx_df <- idx_df %>%
    group_by(Var1, Pathway) %>%
    summarise(pos = sum(value > 0), neg = sum(value < 0))
idx_df$keep <- (idx_df$pos >= 5) | (idx_df$neg >= 5)
radj <- radj[unique(idx_df[idx_df$keep, ]$Var1), ]
rownames(radj) <- gaude_trans[rownames(radj), ]$updated

tmp.text <- gpar(fontfamily="sans", fontsize = 8)
myCols <- rev(brewer.pal(n=5, name="RdBu"))[c(1, 3, 5)]
column_annot <- rep(dataset$DataSets, each = 4)
pdf(paste0("plot/", dir.used, file.prefix, "heatmap_by_Datasets_changed_color.pdf"), width = 8, height = 6)
Heatmap(radj, col = circlize::colorRamp2(c(-0.5, 0, 0.5), myCols), 
        rect_gp = gpar(col = "white"), 
        cluster_columns = FALSE, cluster_rows = TRUE, 
        row_names_side = "left", show_row_dend = FALSE, 
        row_names_max_width = unit(10, "cm"),
        row_names_gp = tmp.text, column_names_gp = tmp.text, 
        row_title_gp = tmp.text, column_title_gp = tmp.text, 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        top_annotation = columnAnnotation(Tumor = column_annot, col = list(Tumor = dataset_color), 
                                          simple_anno_size = unit(2, "mm"), 
                                          annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text)), 
        column_split = column_annot, gap = unit(1, "mm"), 
        column_labels = sub("_.+", "", colnames(radj)), column_names_side = "top",
        column_names_max_height = unit(6, "cm"),
        border = TRUE)
dev.off()

# each pathway need 1/9 inch height, each signature need 1.5cm width
pdf(paste0("plot/", dir.used, file.prefix, "heatmap_by_Sig_changed_color.pdf"), width = 7, height = 5)
Heatmap(radj, col = circlize::colorRamp2(c(-0.5, 0, 0.5), myCols), 
        rect_gp = gpar(col = "white"), 
        cluster_columns = FALSE, cluster_rows = TRUE, 
        row_names_side = "left", show_row_dend = FALSE, 
        row_names_max_width = unit(10, "cm"),
        row_names_gp = tmp.text, column_names_gp = tmp.text, 
        column_title_rot = 0, height = unit(nrow(radj)/9, "inches"), width = unit(1.5*ncol(radj)/10, "cm"), 
        row_title_gp = tmp.text, column_title_gp = tmp.text, 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        top_annotation = HeatmapAnnotation(Tumor = column_annot, col = list(Tumor = dataset_color), 
                                           simple_anno_size = unit(2, "mm"), 
                                           annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
                                           which = "column"), 
        column_split = sub("_.+", "", colnames(radj)), gap = unit(1, "mm"), 
        show_column_names = FALSE, column_names_side = "top", 
        column_names_max_height = unit(6, "cm"),
        border = TRUE)
dev.off()

