# ====================================================================
# Author: Zhou Zhe
# Function: subtype compare of BRCA and COAD
# Version: 1.0
# Date: Aug 18, 2022
# ====================================================================

rm(list=ls())
library(Seurat)
library(scMetab)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(ggpubr)

dir.used <- "3_epi/InterSample/"
dir.create(file.path("res", dir.used), recursive = T)
dir.create(file.path("plot", dir.used), recursive = T)
# file.prefix <- "patient_celltype"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
seed.use = 2021


# I. Load data ====================================================================================
pooled.merged <- readRDS(file = "data/rds/pseudo_sam_celltype.rds")
pooled.merged <- subset(pooled.merged, celltype == "Malignant" & n >= 20)
metab.gsva <- readRDS("data/rds/pseudo_sam_celltype_metab_gsva.rds")
metab.gsva <- metab.gsva[, colnames(pooled.merged)]
pooled.merged[["GSVA"]] <- CreateAssayObject(counts = metab.gsva)

sam_clin <- read.delim("data/Pooled_clinical_info.tsv", na.strings = c("NA", ""))
samCNV <- read.delim(file = "res/3_epi/InterSample/sam_CNV_ITH.tsv")
sam_clin <- merge(sam_clin, samCNV, by = "samID", all = T)
sam_clin <- sam_clin[match(pooled.merged$samID, sam_clin$samID), ]
rownames(sam_clin) <- paste(sam_clin$dataset, sam_clin$samID, "Malignant", sep = "_")


# II. BRCA ========================================================================================
brca_clin <- subset(sam_clin, cancerType == "BRCA" & subtype != "HER2_BC")
brca_gsva <- data.frame(t(metab.gsva[, rownames(brca_clin)]), check.names = F)

res <- plyr::ldply(brca_gsva, function(x) {
    tres = t.test(x[brca_clin$subtype == "ER_BC"], x[brca_clin$subtype == "TNBC"])
    data.frame(statistic = tres$statistic, pvalue = tres$p.value)
})
res$padj <- p.adjust(res$pvalue, method = "BH")

tumorType_color <- celltype_color[c(1, 5)]
names(tumorType_color) <- c("TNBC", "ERBC")

# Horizontal
brca_plot <- subset(res, padj < 0.05)
brca_plot <- brca_plot[order(brca_plot$statistic), ]
brca_plot$.id <- gaude_trans[brca_plot$.id, ]$updated
brca_plot$.id <- factor(brca_plot$.id, levels = brca_plot$.id)
brca_plot$type <- ifelse(brca_plot$statistic > 0, "ERBC", "TNBC")

tmp.text <- element_text(family="sans", size=8)
p1 <- ggplot(brca_plot, mapping = aes(x = .id, y = statistic, fill = statistic < 0)) + 
    geom_bar(stat = "identity", width = 0.6) +
    # coord_flip() +
    # scale_fill_manual(values = ER_TNBC_color) + 
    ggtitle("ER_BC vs TNBC") + 
    xlab("Pathway") + ylab("t value") + 
    theme_classic() + 
    theme(text = tmp.text, axis.text = tmp.text, legend.position = "none", plot.title = tmp.text,
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.line=element_line(size = .3, colour="black"), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 2), "cm"))
ggsave(p1, filename = paste0("plot/3_epi/InterSample/BRCA_subtype_hist_horizontal.pdf"), width = 6, height = 5)

# vertical
tmp.text <- element_text(family="sans", size=8)
p2 <- ggplot(brca_plot, mapping = aes(x = .id, y = statistic, fill = type)) + 
    geom_bar(stat = "identity", width = 0.6) +
    coord_flip() +
    # scale_fill_manual(values = brewer.pal(2, "RdYlBu")) + 
    ggtitle("ER_BC vs TNBC") + 
    xlab("Pathway") + ylab("t value") + 
    scale_fill_manual(values = tumorType_color) +
    ylim(-10, 5) +
    theme_classic() + 
    theme(text = tmp.text, axis.text = tmp.text, legend.position = "none", plot.title = tmp.text,
          # axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.line=element_line(size = .3, colour="black"), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
ggsave(p2, filename = paste0("plot/3_epi/InterSample/BRCA_subtype_hist_vertival.pdf"), width = 5, height = 4)

# boxplot of specific pathway
tumorType_color <- celltype_color[c(1, 5)]
names(tumorType_color) <- c("TNBC", "ER_BC")
Imp_path = c("Glycolysis and Gluconeogenesis", "Carnitine shuttle", "Ketone Bodies Metabolism")
pdf(file = paste0("plot/3_epi/InterSample/BRCA_subtype_boxplot_Imp.pdf"), width = 4, height = 2)
plot_data <- cbind(brca_gsva, brca_clin)
for (path in Imp_path) {
    p = ggplot(plot_data, mapping = aes(x = subtype, y = .data[[path]]), color = subtype) + 
        geom_boxplot(aes(color = subtype), width = 0.3) + 
        geom_jitter(position = position_jitter(0.2), aes(color = subtype)) +  
        scale_color_manual(values = tumorType_color) + 
        theme_classic() + 
        facet_wrap(~dataset) +
        stat_compare_means(method = "t.test", size = 2.88, 
                           label.y.npc = 0.9, 
                           aes(label = paste0("P = ", after_stat(p.format)))) + 
        theme(text = tmp.text, axis.text = tmp.text, legend.position = "none", plot.title = tmp.text,
              strip.text = tmp.text, 
              axis.line=element_line(size = .3, colour="black"),
              plot.margin = unit(c(0.2, 0.5, 0.5, 0.5), "cm"))
    print(p)
}
dev.off()



# II. COAD ========================================================================================
coad_clin <- subset(sam_clin, cancerType == "COAD")
coad_gsva <- data.frame(t(metab.gsva[, rownames(coad_clin)]), check.names = F)

res <- plyr::ldply(coad_gsva, function(x) {
    tres = t.test(x[coad_clin$subtype == "MMRd"], x[coad_clin$subtype == "MMRp"])
    data.frame(statistic = tres$statistic, pvalue = tres$p.value)
})
res$padj <- p.adjust(res$pvalue, method = "BH")

# Horizontal
coad_plot <- subset(res, padj < 0.05)
coad_plot <- coad_plot[order(coad_plot$statistic), ]
coad_plot$.id <- factor(coad_plot$.id, levels = coad_plot$.id)
coad_plot$type <- ifelse(coad_plot$statistic > 0, "MMRd", "MMRp")

tumorType_color <- celltype_color[c(1, 5)]
names(tumorType_color) <- c("MMRp", "MMRd")

# vertical
tmp.text <- element_text(family="sans", size=8)
p2 <- ggplot(coad_plot, mapping = aes(x = .id, y = statistic, fill = type)) + 
    geom_bar(stat = "identity", width = 0.6) +
    coord_flip() +
    # scale_fill_manual(values = brewer.pal(2, "RdYlBu")) + 
    ggtitle("MMRd vs MMRp") + 
    xlab("Pathway") + ylab("t value") + 
    scale_fill_manual(values = tumorType_color) +
    theme_classic() + 
    theme(text = tmp.text, axis.text = tmp.text, legend.position = "none", plot.title = tmp.text,
          # axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.line=element_line(size = .3, colour="black"), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
ggsave(p2, filename = paste0("plot/3_epi/InterSample/COAD_subtype_hist_vertival.pdf"), width = 4, height = 2.5)

# boxplot of specific pathway
Imp_path = c("Glycolysis and Gluconeogenesis", "Fatty Acid Biosynthesis", "Citric Acid Cycle")
pdf(file = paste0("plot/3_epi/InterSample/COAD_subtype_boxplot_Imp.pdf"), width = 1.9, height = 2)
plot_data <- cbind(coad_gsva, coad_clin)
for (path in Imp_path) {
    p = ggplot(plot_data, mapping = aes(x = subtype, y = .data[[path]]), color = subtype) + 
        geom_boxplot(aes(color = subtype), width = 0.3) + 
        geom_jitter(position = position_jitter(0.2), aes(color = subtype)) +  
        scale_color_manual(values = tumorType_color) + 
        theme_classic() +
        stat_compare_means(method = "t.test", size = 2.88, 
                           label.y.npc = 0.9, 
                           aes(label = paste0("P = ", after_stat(p.format)))) + 
        theme(text = tmp.text, axis.text = tmp.text, legend.position = "none", plot.title = tmp.text,
              strip.text = tmp.text, 
              axis.line=element_line(size = .3, colour="black"))
    print(p)
}
dev.off()
