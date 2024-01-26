# ====================================================================
# Author: Zhou Zhe
# Function: ROGUE indexes are used to measure homogeneity of cell clusters, this script calculated metabolic
# genes based ROGUE index under two levels:
#     1. sample-wise: calculated for each dataset and each cell type, show heterogeneity between samples
#     2. dataset-wise: calculated for each cell type, show heterogeneity between datasets
# Version: 1.0
# Date: Feb 28, 2023
# ====================================================================

rm(list=ls())
library(ROGUE)
library(Seurat)
library(scMetab)
library(ggplot2)


dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))

# sample level ROGUE index across dataset
seu_obj <- readRDS("data/rds/all_dataset_harmony_metab_features.rds")
metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu_obj))
expr_all <- GetAssayData(seu_obj, slot = "counts")[metab.gene, ]

expr <- matr.filter(as.matrix(expr_all))
meta <- seu_obj[[c("celltype", "patientID", "dataset")]][colnames(expr), ]

res1 <- rogue(expr, labels = meta$celltype, samples = meta$dataset, platform = "UMI", 
             span = 0.6)
rogue.boxplot(res1)
write.csv(res1, file = "res/1_clustering/MetabGene_based_sample_level_ROGUE_index.csv", 
          quote = F)

res1 <- read.csv("res/1_clustering/MetabGene_based_sample_level_ROGUE_index.csv", 
                 row.names = 1, check.names = F)
res1$dataset <- rownames(res1)
tmp_plot <- reshape2::melt(res1, id.vars = "dataset", variable.name = "celltype", value.name = "value")
tmp_plot$celltype <- factor(tmp_plot$celltype, levels = names(celltype_color))

tmp.text <- element_text(family="sans", size=12)
ggplot(tmp_plot, mapping = aes(x = celltype, y = value, group = dataset, colour = dataset)) + 
    geom_point() + 
    geom_line(linewidth = 1) + 
    scale_color_manual(values = dataset_color) + 
    theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text, 
          axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "plot/1_clustering/sample_ROGUE_lineplot.pdf", width = 3, height = 3)


pvalue <- tmp_plot %>% rstatix::anova_test(value ~ celltype)
ggplot(tmp_plot, mapping = aes(x = celltype, y = value, colour = celltype)) + 
    geom_boxplot(width = 0.6) + geom_jitter() + 
    scale_color_manual(values = celltype_color) + 
    theme_classic() + 
    ggtitle("p = 4.5e-11") + 
    theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text, 
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(filename = "plot/1_clustering/sample_ROGUE_boxplot.pdf", width = 7, height = 4)

# dataset level ROGUE index
meta$platform = "UMI"
res2 <- rogue(expr, labels = meta$celltype, samples = meta$platform, platform = "UMI", span = 0.6)
write.csv(res2, file = "res/1_clustering/MetabGene_based_dataset_level_ROGUE_index.csv", 
          quote = F)

tmp.text <- element_text(family="sans", size=6)
tmp_plot2 <- data.frame(t(res2[, c(1,3, 5:8)]), check.names = F)
tmp_plot2$celltype <- factor(rownames(tmp_plot2), 
                             levels = c(intersect(names(celltype_color), rownames(tmp_plot2))))
tmp_plot2$rogue = signif(tmp_plot2$UMI, digits = 2)
ggplot(tmp_plot2, aes(x = celltype, y = UMI)) + 
    geom_segment(aes(x = celltype, xend = celltype, y = 0.9, yend = UMI), 
                 linewidth = 1.8, color = "grey") + 
    geom_point(aes(color = celltype, fill = celltype), size = 4, pch = 21) + 
    geom_text(mapping = aes(x = celltype, y = UMI, label = rogue), size = 2) + 
    scale_fill_manual(values = celltype_color) + 
    scale_color_manual(values = celltype_color) + 
    ylab("ROGUE index") + 
    theme_classic() + 
    theme(text = tmp.text, legend.text = tmp.text, axis.text = tmp.text, 
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(filename = "plot/1_clustering/dataset_ROGUE.pdf", width = 1.2, height = 1.8)
