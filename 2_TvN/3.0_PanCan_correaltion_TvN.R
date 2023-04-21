# ==================================================================================================
# Author: Zhou Zhe
# Program: compare cell-to-cell pca correlation of tumor to normal
#          到底要比较样本间差异还是细胞间差异？再思考一下！
#          Datasets: Pan Cancer
#          TorN: Tumor and Normal
#          Celltype: All celltype
# Version: 1.0
# Date: Aug 16, 2022
# ==================================================================================================

# init -------
rm(list=ls())

library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggsci)
library(scMetab)
library(RColorBrewer)
library(patchwork)
library(parallel)

n_cores <- 8
seed.use = 2021
dir.used <- "2_TvN/TvN_correlation/"
dir.create(file.path("plot/", dir.used))
dir.create(file.path("res/", dir.used))

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
# cl <- makeCluster(10)


plot_list = list()

for (i in c(1, 4, 7:10)) {
    tumor <- dataset$DataSets[i]
    print(tumor)
    
    # Load data ------
    seu_obj <- readRDS(paste0("data/rds/MetabGene/", tumor, ".TpN.unintegrated.rds"))
    seu_obj$celltype_x_TorN <- paste(seu_obj$celltype, seu_obj$TorN, sep = "_")
    
    # Construct compare data.frame ------
    celltypes <- c("Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")
    compare.df <- data.frame(Tumor = paste0(celltypes, "_T"), Normal = paste0(celltypes, "_N"))
    compare.df <- rbind(data.frame(Tumor = "Malignant_T", Normal = "Epithelial_N"), compare.df)
    rownames(compare.df) <- c("Malignant", celltypes)
    
    cell_pairs <- sample_cell_pairs(seu_obj[, seu_obj$celltype_x_TorN %in% unlist(compare.df)], 
                                    group.by = "patientID", split.by = "celltype_x_TorN",
                                    n = 500, seed.use = seed.use)
    
    metab.pca <- t(Embeddings(seu_obj, reduction = "pca"))[1:30, ]
    pca_cor <- cell_pairs
    # clusterExport(cl, varlist = list("metab.pca"))
    system.time(
        pca_cor$cor <- apply(cell_pairs, 1, function(x) {
            cor(metab.pca[, x[1]], metab.pca[, x[2]], method = "spearman")
        })
    )
    # stopCluster(cl)
    
    
    pca_cor$celltype <- sub("_.+", "", pca_cor$split.by)
    pca_cor$TorN <- sub(".+_", "", pca_cor$split.by)
    pca_cor$celltype <- ifelse(pca_cor$celltype %in% c("Malignant", "Epithelial"),
                               yes = "Malignant/Epithelial", no = pca_cor$celltype)
    
    tmp.text <- element_text(family="sans", size=8)
    p1 <- ggplot(data = pca_cor[pca_cor$pair == "inter", ], mapping = aes(x = celltype, y = cor, fill = TorN)) +
        geom_split_violin(trim = FALSE, color = "white", width = 1, scale = "width") +
        geom_boxplot(width = 0.15, position = position_dodge(0.3),
                     coef = 0, outlier.shape = NA, lwd = 0.2) +
        ylim(-0.5, 1) + xlab("Cell Type") + ylab("Spearman correlation") + ggtitle(tumor) +
        theme(text=tmp.text, axis.text = tmp.text, legend.text = tmp.text, plot.title = tmp.text, 
              axis.line=element_line(size = .3, colour="black"),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = TvN_color)
    # tmpfile <- paste0("plot/", dir.used, file.prefix, "intersample_cor_distribution_pca.pdf")
    # ggsave(p, file = tmpfile, width = 4, height = 2)
    
    
    p2 <- ggplot(data = pca_cor[pca_cor$pair == "intra", ], mapping = aes(x = celltype, y = cor, fill = TorN)) +
        geom_split_violin(trim = FALSE, color = "white", width = 1, scale = "width") +
        geom_boxplot(width = 0.15, position = position_dodge(0.3),
                     coef = 0, outlier.shape = NA, lwd = 0.2) +
        ylim(-0.5, 1.5) + xlab("Cell Type") + ylab("Spearman correlation") + ggtitle(tumor) + 
        theme(text=tmp.text, axis.text = tmp.text, legend.text = tmp.text, plot.title = tmp.text, 
              axis.line=element_line(size = .3, colour="black"),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = TvN_color)
    # tmpfile <- paste0("plot/", dir.used, file.prefix, "intrasample_cor_distribution_pca.pdf")
    # ggsave(p, file = tmpfile, width = 4, height = 2.5)
    
    tmpfile <- paste0("res/", dir.used, tumor, "_cor_distribution_pca.tsv")
    write.table(pca_cor, file = tmpfile, quote = F, sep = "\t", row.names = F)
    
    plot_list[[tumor]] <- list(inter = p1, intra = p2)
}



tmpfile <- paste0("plot/", dir.used, file.prefix, "TvN_cor_distribution_pca.pdf")
ggsave(p1 + p2, file = tmpfile, width = 6, height = 2.5)

###### ====== Load Data ====== ######

render_plot <- function(p) {
    p[[1]] <- p[[1]] + xlab(NULL) + theme(axis.text.x = element_blank(), legend.position = "none")
    p[[2]] <- p[[2]] + ylab(NULL) + xlab(NULL) + 
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")
    p[[3]] <- p[[3]] + ylab(NULL) + xlab(NULL) +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    p[[4]] <- p[[4]] + xlab(NULL) + theme(legend.position = "none")
    p[[5]] <- p[[5]] + ylab(NULL) + theme(axis.text.y = element_blank(), legend.position = "none")
    p[[6]] <- p[[6]] + xlab(NULL) + ylab(NULL) + theme(axis.text.y = element_blank())
    return(wrap_plots(p))
}

ggsave(render_plot(lapply(plot_list, function(x) x$inter)), 
       file = paste0("plot/", dir.used, "PanCan_TvN_Intersample_cor_distribution_pca.pdf"), 
       width = 6, height = 3)
ggsave(render_plot(lapply(plot_list, function(x) x$intra)), 
       file = paste0("plot/", dir.used, "PanCan_TvN_Intrasample_cor_distribution_pca.pdf"), 
       width = 6, height = 3)

##### ====== pca spearman correlation ====== ######






