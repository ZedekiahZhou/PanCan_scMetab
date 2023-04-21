# ====================================================================
# Author: Zhou Zhe
# Function: Find inter-sample variable metabolic pathways for malignant 
#           Use pooled data from 6_patient_clustering/1.0_patient_pooled.R
# Version: 1.0
# Date: Aug 17, 2022
# ====================================================================

rm(list=ls())
library(Seurat)
library(scMetab)
library(ggplot2)

dir.used <- "3_epi/InterSample/"
dir.create(file.path("res", dir.used), recursive = T)
dir.create(file.path("plot", dir.used), recursive = T)
# file.prefix <- "patient_celltype"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
seed.use = 2021

# I. Inter sample variable pathway ================================================================
## 1. Load poolted data ----
pooled.merged <- readRDS(file = "data/rds/pseudo_sam_celltype.rds")
pooled.merged <- subset(pooled.merged, celltype == "Malignant" & n >= 20)
metab.gsva <- readRDS("data/rds/pseudo_sam_celltype_metab_gsva.rds")
metab.gsva <- metab.gsva[, colnames(pooled.merged)]
pooled.merged[["GSVA"]] <- CreateAssayObject(counts = metab.gsva)

DefaultAssay(pooled.merged) <- "GSVA"
pooled.list <- SplitObject(pooled.merged, split.by = "tumor")
pooled.list <- lapply(pooled.list, FindVariableFeatures)

## 2. top 10 path ------
intersam_feature <- lapply(1:10, function(i) {
    tumor <- dataset$DataSets[i]
    tmp <- pooled.list[[tumor]]@assays$GSVA@meta.features
    tmp$pathway <- rownames(tmp)
    tmp <- tmp[order(tmp$vst.variance.standardized, decreasing = TRUE)[1:10], ]
    tmp$dataset <- tumor
    tmp
})
intersam_feature <- do.call(rbind, intersam_feature)

## 3. circle hist plot ------
intersam_feature$id <- 1:nrow(intersam_feature)
angle <- 90 - 360 * (intersam_feature$id - 0.5) / nrow(intersam_feature)
intersam_feature$hjust <- ifelse(angle < -90, 1, 0)
intersam_feature$angle <- ifelse(angle < -90, angle + 180, angle)
intersam_feature$pathway <- gaude_trans[intersam_feature$pathway, ]$updated

tmpfile <- paste0("plot/", dir.used, "PanCan_inter_sample_path_change_color.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 7, height = 8)
p <- ggplot(intersam_feature, aes(x = id, y = vst.variance.standardized, fill = dataset)) +
    geom_bar(stat = "identity", alpha = 1) +
    ylim(-2, max(intersam_feature$vst.variance.standardized)+30) +
    theme_minimal() +
    scale_fill_manual(values = dataset_color) +
    theme(
        legend.position = "top",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        text = tmp.text
    ) +
    coord_polar() +
    geom_text(data = intersam_feature, aes(x = id, y = max(vst.variance.standardized)+1, label = pathway, hjust = hjust),
              color = "black", size = 2.82, angle = intersam_feature$angle, inherit.aes = FALSE)
p
dev.off()

## 4. stacked circile hist --------
label_data <- intersam_feature %>%
    group_by(pathway) %>%
    summarise(n = n(), stats = sum(vst.variance.standardized))
label_data <- label_data[order(label_data$n, label_data$stats, decreasing = TRUE), ]
intersam_feature$pathway <- factor(intersam_feature$pathway, levels = label_data$pathway)
label_data$id <- 1:nrow(label_data)
angle <- 90 - 360 * (label_data$id - 0.5) / nrow(label_data)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)

tmpfile <- paste0("plot/", dir.used, "PanCan_inter_sample_stacked_change_color.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 6, height = 8)
p <- ggplot(intersam_feature, aes(x = pathway, y = vst.variance.standardized, fill = dataset)) + 
    geom_bar(stat = "identity", alpha = 1) + 
    ylim(-10, max(label_data$stats) + 31) + 
    theme_minimal() +
    scale_fill_manual(values = dataset_color) + 
    theme(
        legend.position = "bottom", 
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        text = tmp.text
    ) +
    coord_polar() + 
    geom_text(data = label_data, aes(x = id, y = stats + 2, label = pathway, hjust = hjust), 
              color = "black", size = 2.82, angle = label_data$angle, inherit.aes = FALSE)
p
dev.off()



