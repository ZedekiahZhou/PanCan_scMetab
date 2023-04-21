# ==================================================================================================
# Author: Zhou Zhe
# Program: Compare similarity of cells (inter-patients vs intra-patients)
#          Methods: PCA cor based
#          TorN: Tumor only
#          Celltype: All celltype included
# Notes: metabolic genes correlation will be impact by nFeature_metab
#        use PCA correlation!!!
# Version: 1.0
# Date: Aug 14, 2022
# ==================================================================================================

# init -------
rm(list=ls())

library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(pheatmap)
library(scMetab)
library(patchwork)

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

n_cores <- 8
seed.use = 2021
dir.used <- "1_clustering/correlation/"
dir.create(file.path("plot/", dir.used))
dir.create(file.path("res/", dir.used))

for (i in c(1, 4, 7:10)) {
    tumor <- dataset$DataSets[i]
    print(tumor)
    
    seu_obj <- readRDS(paste0("data/rds/MetabGene/", tumor, ".TpN.unintegrated.rds"))
    seu_obj <- subset(seu_obj, TorN == "T")
    
    cell_pairs <- sample_cell_pairs(seu_obj, group.by = "patientID", split.by = "celltype",
                                    n = 500, seed.use = seed.use)
    
    metab.pca <- t(Embeddings(seu_obj, reduction = "pca"))[1:30, ]
    pca_cor <- cell_pairs
    system.time(
        pca_cor$cor <- apply(cell_pairs, 1, function(x) {
            cor(metab.pca[, x[1]], metab.pca[, x[2]], method = "spearman")
        })
    )
    
    
    tmpfile <- paste0("res/", dir.used, tumor, "_cor_distribution_pca.tsv")
    write.table(pca_cor, file = tmpfile, quote = F, sep = "\t", row.names = F)
    
    plot_list[[tumor]] = p
}


used_celltypes <- c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")
plot_list <- lapply(c(1, 4, 7:10), function(i) {
    tumor <- dataset$DataSets[i]
    print(tumor)
    
    pca_cor <- read.delim(paste0("res/", dir.used, tumor, "_cor_distribution_pca.tsv"))
    pca_cor <- subset(pca_cor, split.by %in% used_celltypes)
    
    tmp.text <- element_text(family="sans", size=8)
    p <- ggplot(data = pca_cor, mapping = aes(x = split.by, y = cor, fill = pair)) +
        geom_split_violin(trim = FALSE, color = "white", width = 1) +
        geom_boxplot(width = 0.15, position = position_dodge(0.3),
                     coef = 0, outlier.shape = NA, lwd = 0.2) +
        ylim(-0.4, 1) + xlab("Cell Type") + ylab("Spearman correlation") +
        ggtitle(tumor) + 
        theme(text=tmp.text, 
              axis.text = tmp.text, 
              legend.text = tmp.text, 
              plot.title = tmp.text, 
              axis.line=element_line(size = .3, colour="black"),
              axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
})

plot_list[[1]] <- plot_list[[1]] + xlab(NULL) + theme(axis.text.x = element_blank(), legend.position = "none")
plot_list[[2]] <- plot_list[[2]] + ylab(NULL) + xlab(NULL) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")
plot_list[[3]] <- plot_list[[3]] + ylab(NULL) + xlab(NULL) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_list[[4]] <- plot_list[[4]] + theme(legend.position = "none")
plot_list[[5]] <- plot_list[[5]] + ylab(NULL) + theme(axis.text.y = element_blank(), legend.position = "none")
plot_list[[6]] <- plot_list[[6]] + ylab(NULL) + theme(axis.text.y = element_blank())

wrap_plots(plot_list)
tmpfile <- paste0("plot/", dir.used, "PanCan_cor_distribution_pca.pdf")
ggsave(wrap_plots(plot_list), file = tmpfile, width = 6, height = 3)


