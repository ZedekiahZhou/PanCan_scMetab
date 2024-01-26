# ====================================================================
# Author: Zhou Zhe
# Function: Identify Important TF for Malignant cells
# Version: 1.0
# Date: Oct 30, 2022
# ====================================================================

# 0. INIT ============================================================
rm(list=ls())
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggsci)
library(ggvenn)
library(ggrepel)

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
main_dataset <- dataset$DataSets[c(1, 4, 7:10)]

dir.used <- "3_epi/Identify_Imp_TF/"
dir.create(file.path("plot", dir.used), recursive = T)

# # I. Filter core targets ==========================================================================
# #### 首先筛选关键通路共有的targets
# celltype_Metab_markers <- read.delim("res/1_clustering/PAS_markers/Celltype_Metab_Markers_All_Datasets.tsv")
# celltype_Metab_markers <- subset(celltype_Metab_markers, cluster == "Malignant")
# Malig_core <- lapply(unique(celltype_Metab_markers$gene), function(pathway) {
#     genes = celltype_Metab_markers[celltype_Metab_markers$gene == pathway, c("tumor", "core")]
#     genes = separate_rows(genes, core, sep = "/")
#     genes = genes %>%
#         group_by(core) %>%
#         summarise(freq = n(), tumor = paste(tumor, collapse = "/"))
#     genes = genes[order(genes$freq, decreasing = T), ]
# })
# names(Malig_core) <- unique(celltype_Metab_markers$gene)
# 
# TvN_Metab_markers <- read.delim("res/2_TvN/DE/TvN_Metab_Markers_All_Datasets.tsv")
# TvN_Metab_markers <- subset(TvN_Metab_markers, cluster == "Malignant")
# TvN_core <- lapply(unique(TvN_Metab_markers$pathway), function(pathway) {
#     genes = TvN_Metab_markers[TvN_Metab_markers$pathway == pathway, c("tumor", "core")]
#     genes = separate_rows(genes, core, sep = "/")
#     genes = genes %>%
#         group_by(core) %>%
#         summarise(freq = n(), tumor = paste(tumor, collapse = "/"))
#     genes = genes[order(genes$freq, decreasing = T), ]
# })
# names(TvN_core) <- unique(TvN_Metab_markers$pathway)
# 
# # 筛选得到的糖酵解关键targets
# intersect(Malig_core$`Glycolysis and Gluconeogenesis`$core, TvN_core$`Glycolysis and Gluconeogenesis`$core)
# core_targets <- c("ALDOA", "GAPDH", "GPI", "LDHA", "LDHB", "PFKP", "PFKL", "PGK1", "SLC2A1", "ENO1", "ENO2", 
#                   "PKM", "TPI1", "PGAM1", "SLC16A3")
# 
# 
# 
# # II. Filter TF by DE ===================================================================================
# 
# ## 调整筛选思路，由于MMP的重点是肿瘤内异质性，所以MMP的TF通过肿瘤相较于其他细胞高表达来看并不合适！！
# 
# 
# ## Common TF markers for Malignant
# celltype_TF_markers <- read.delim("res/01_clustering/PAS_markers/Celltype_TF_Markers_All_Datasets.tsv")
# celltype_TF_markers <- subset(celltype_TF_markers, cluster == "Malignant")
# celltype_TF_pos <- subset(celltype_TF_markers, avg_log2FC > 0)
# Malig_top_TF <- sort(table(celltype_TF_pos$gene), decreasing = T)
# Malig_top_TF <- as.data.frame(Malig_top_TF[Malig_top_TF >=4])
# Malig_top_TF$Var1 <- paste0(Malig_top_TF$Var1, "_reg")
# Malig_top_TF$reg <- factor(Malig_top_TF$Var1, levels = rev(Malig_top_TF$Var1))
# 
# tmp.text <- element_text(family="sans", size=8)
# p1 <- ggplot(as.data.frame(Malig_top_TF[1:10, ]), mapping = aes(x = reg, y = Freq, fill = "A")) + 
#     geom_bar(stat = "identity", width = 0.5) +
#     coord_flip() +
#     scale_fill_manual(values = brewer.pal(11, "RdYlBu")[2]) + 
#     ggtitle("Malignant Top Regulons") + 
#     xlab("Regulons") + ylab("Number of datasets") +
#     theme_classic() + 
#     theme(text = tmp.text, axis.text = tmp.text, legend.position = "none", plot.title = tmp.text,
#           axis.line=element_line(size = .3, colour="black"))
# 
# ## Common TF markers for TvN
# TvN_TF_markers <- read.delim("res/03_TvN/DE/TvN_TF_Markers_All_Datasets.tsv")
# TvN_TF_markers <- subset(TvN_TF_markers, cluster == "Malignant")
# TvN_TF_pos <- subset(TvN_TF_markers, avg_log2FC > 0)
# TvN_top_TF <- sort(table(TvN_TF_pos$pathway), decreasing = T)
# TvN_top_TF <- as.data.frame(TvN_top_TF[TvN_top_TF >= 4])
# TvN_top_TF$Var1 <- paste0(TvN_top_TF$Var1, "_reg")
# TvN_top_TF$reg <- factor(TvN_top_TF$Var1, levels = rev(TvN_top_TF$Var1))
# 
# tmp.text <- element_text(family="sans", size=8)
# p2 <- ggplot(as.data.frame(TvN_top_TF[1:10, ]), mapping = aes(x = reg, y = Freq, fill = "A")) + 
#     geom_bar(stat = "identity", width = 0.5) +
#     coord_flip() +
#     scale_fill_manual(values = brewer.pal(11, "RdYlBu")[2]) + theme_classic() +
#     ggtitle("Tumor vs. Normal Top Regulons") + 
#     xlab("Regulons") + ylab("Number of datasets") +
#     theme(text = tmp.text, axis.text = tmp.text, legend.position = "none", plot.title = tmp.text,
#           axis.line=element_line(size = .3, colour="black"))
# 
# ggsave(p1 + p2, filename = paste0("plot/", dir.used, "Malignant_topTFs.pdf"), width = 6, height = 2)
# 
# p3 <- ggvenn(list(Malig_top_TF = as.character(Malig_top_TF$Var1), TvN_top_TF = as.character(TvN_top_TF$Var1)), 
#        fill_color = brewer.pal(9, "Set1")[1:2], show_percentage = FALSE, fill_alpha = 0.5, 
#        set_name_size = 2.82, text_size = 2.82)
# ggsave(p3, filename = paste0("plot/", dir.used, "TF_venn.pdf"), width = 2, height = 2)
# 
# 
# # comman TF
# con_TF <- sub("_reg", "", intersect(Malig_top_TF$Var1, TvN_top_TF$Var1))
# 
# # TFs regulate glucolysis genes 
# TvN_TF_targets <- TvN_TF_markers[TvN_TF_markers$pathway %in% con_TF, c("pathway", "tumor", "core")]
# TvN_TF_targets <- separate_rows(TvN_TF_targets, core, sep = "/")
# TvN_TF_targets <- TvN_TF_targets[TvN_TF_targets$core %in% core_targets, ]
# target_counts <- as.data.frame.array(table(TvN_TF_targets$pathway, TvN_TF_targets$core))
# target_counts$n = rowSums(target_counts > 0)
# 
# #### Plot conTF regulons for TvN 
# TF_reg <- subset(TvN_TF_markers, pathway %in% con_TF)
# TF_reg <- dcast(TF_reg, pathway ~ tumor, value.var = "avg_log2FC") %>%
#     column_to_rownames(var = "pathway")
# TF_reg[is.na(TF_reg)] <- 0
# 
# tmp.text <- gpar(fontfamily="sans", fontsize = 8)
# tmpfile <- paste0("plot/", dir.used, "Malignant_conTF_reg_TvN.pdf")
# pdf(tmpfile, width = 5, height = 3)
# Heatmap(TF_reg, circlize::colorRamp2(seq(0.1, -0.1, -0.05), brewer.pal(5, "RdBu")), 
#         rect_gp = gpar(col = "white"), 
#         row_labels = paste0(rownames(TF_reg), "_reg"),
#         row_names_gp = tmp.text, column_names_gp = tmp.text, 
#         column_names_rot = 45, 
#         row_title_gp = tmp.text, column_title_gp = tmp.text, 
#         heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
#         name = "log2FC", 
#         cluster_rows = FALSE, cluster_columns = FALSE)
# dev.off()
# 
# #### Plot conTF expression for TvN
# TvN_gene_markers <- read.delim(file = paste0("res/2_TvN/DE/TvN_Gene_Markers_All_Datasets.tsv"))
# TF_expr <- subset(TvN_gene_markers, gene %in% con_TF & cluster == "Malignant")
# TF_expr <- dcast(TF_expr, gene~tumor, value.var = "avg_log2FC") %>%
#     column_to_rownames(var = "gene")
# TF_expr[is.na(TF_expr)] <- 0
# 
# tmp.text <- gpar(fontfamily="sans", fontsize = 8)
# tmpfile <- paste0("plot/", dir.used, "Malignant_conTF_expr_TvN.pdf")
# pdf(tmpfile, width = 5, height = 3)
# Heatmap(TF_expr, circlize::colorRamp2(seq(1, -1, -0.5), brewer.pal(5, "RdBu")), 
#         rect_gp = gpar(col = "white"), 
#         row_names_gp = tmp.text, column_names_gp = tmp.text, 
#         column_names_rot = 45, 
#         row_title_gp = tmp.text, column_title_gp = tmp.text, 
#         heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
#         name = "log2FC", 
#         cluster_rows = FALSE, cluster_columns = FALSE)
# dev.off()


# III. Filter TF by survival ====================================================================
MMP_TF <- read.delim("res/02_NMF/TF_of_MMPs.tsv")
path.survival <- read.csv(file = "res/08_patient_clustering/PathSurv/path.survival.csv")


## TF expr and survival 
TF.survival = subset(path.survival, dataset %in% c("BRCA", "LUAD", "COAD", "PAAD", "STAD") & 
               clu_method %in% rownames(MMP_TF))
TF.survival$Sig = ifelse(TF.survival$Pval >= 0.05, "No_effect", 
                         ifelse(TF.survival$HR > 1, "Better", "Worse"))
TF.survival$label = ifelse(TF.survival$Sig != "No_effect", TF.survival$clu_method, "")
Imp_TF <- c("CEBPG", "E2F2", "EZH2", "FOSL1", "KLF3", "MYBL2")  # TF survival significant in at least 2 cancers
TF.survival$label2 = ifelse(TF.survival$clu_method %in% Imp_TF & TF.survival$Sig != "No_effect", TF.survival$clu_method, "")

tmpcolor <- brewer.pal(9, "Set1")[c(1,2,11)]
names(tmpcolor) = c("Worse", "Better", "No_effect")
tmp.text <- element_text(family="sans", size=6)
ggplot(TF.survival, aes(x = cancer_type, y = log2(HR), color = Sig, label = label2)) + 
    geom_jitter(shape = 1, size = 1.5/.pt, position = position_jitter(seed = 1, width = 0.1)) + 
    scale_color_manual(values = tmpcolor) + 
    geom_text_repel(max.overlaps = 100, show.legend = F, size = 6/.pt, segment.size = 0.3, 
                    min.segment.length = 0, position = position_jitter(seed = 1, width = 0.1)) +
    coord_flip() + theme_classic() + xlab("Cancer type") + ylab("Log2(Hazard ratio)") + 
    theme(axis.title = tmp.text, axis.text = tmp.text, legend.text = tmp.text, legend.title = tmp.text, 
          legend.position = "top")
ggsave("plot/04_epi/Identify_Imp_TF/TF_survival_HR.pdf", width = 3, height = 2)






# TF.survival <- subset(path.survival, dataset %in% c("BRCA", "LUAD", "COAD", "PAAD", "STAD") & 
#                           clu_method %in% con_TF)
# idx <- TF.survival$Pval > 0.1
# TF.survival$Survival <- ifelse(TF.survival$Pval > 0.1, "No_effect", 
#                              ifelse(TF.survival$HR < 1, "Worse", "Better"))
# 
# # tmpfile <- paste0("plot/", dir.used, "ALL_pooled_dotplot.pdf")

# pdf("plot/2_TvN/Identify_Imp_TF/TF_expr_survival_dotplot.pdf", width = 6, height = 2.5)
# # pdf(tmpfile, width = 12, height = length(levels(TvN_metab_markers$pathway))/9)
# print(ggplot(TF.survival, aes(x = clu_method, y = cancer_type, color = Survival)) + 
#           geom_point(aes(size = -log10(Pval))) + 
#           scale_color_manual(values = c(Worse = tmpcolor[1], Better = tmpcolor[2], No_effect = tmpcolor[3])) + 
#           scale_size(limits = c(0, 4), range = c(0, 4)) + 
#           theme_classic() + 
#           xlab("Pathway") + ylab("Cancer type") + 
#           theme(axis.text=tmp.text, text = tmp.text, 
#                 axis.text.x = element_text(angle = 45, hjust = 1), 
#                 axis.line=element_line(size = .3, colour="black"), 
#                 legend.position = "top", 
#                 plot.margin = margin(0.5, 2, 0.5, 2, unit = "cm")))
# dev.off()
# 
# 
# ## TF reg and survival 
# reg.survival <- subset(path.survival, dataset %in% c("BRCA", "LUAD", "COAD", "PAAD", "STAD") & 
#                           clu_method %in% paste0(con_TF, "_reg"))
# idx <- reg.survival$Pval > 0.1
# reg.survival$Survival <- ifelse(reg.survival$Pval > 0.1, "No_effect", 
#                                ifelse(reg.survival$HR < 1, "Worse", "Better"))
# 
# # tmpfile <- paste0("plot/", dir.used, "ALL_pooled_dotplot.pdf")
# tmpcolor <- brewer.pal(9, "Set1")[c(1,2,11)]
# pdf("plot/2_TvN/Identify_Imp_TF/TF_reg_survival_dotplot.pdf", width = 6, height = 2.5)
# # pdf(tmpfile, width = 12, height = length(levels(TvN_metab_markers$pathway))/9)
# print(ggplot(reg.survival, aes(x = clu_method, y = cancer_type, color = Survival)) + 
#           geom_point(aes(size = -log10(Pval))) + 
#           scale_color_manual(values = c(Worse = tmpcolor[1], Better = tmpcolor[2], No_effect = tmpcolor[3])) + 
#           scale_size(limits = c(0, 4), range = c(0, 4)) + 
#           theme_classic() + 
#           xlab("Pathway") + ylab("Cancer type") + 
#           theme(axis.text=tmp.text, text = tmp.text, 
#                 axis.text.x = element_text(angle = 45, hjust = 1), 
#                 axis.line=element_line(size = .3, colour="black"), 
#                 legend.position = "top", 
#                 plot.margin = margin(0.5, 2, 0.5, 2, unit = "cm")))
# dev.off()
