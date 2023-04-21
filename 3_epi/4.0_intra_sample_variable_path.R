# ====================================================================
# Author: Zhou Zhe
# Function: Find intra-sample variable metabolic pathways for malignant 
#           Use pooled data from 6_patient_clustering/1.0_patient_pooled.R
# Version: 1.0
# Date: Aug 17, 2022
# ====================================================================

rm(list=ls())
library(Seurat)
library(scMetab)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)

dir.used <- "3_epi/IntraSample/"
dir.create(file.path("res", dir.used), recursive = T)
dir.create(file.path("plot", dir.used), recursive = T)
# file.prefix <- "patient_celltype"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
seed.use = 2021


# I. calculate CNV score and ITH ==================================================================
intrasam_feature = list()
for (i in 1:10) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    seu_obj <- readRDS(paste0("data/rds/PAS//", tumor, ".TpN.PAS.unintegrated.rds"))
    DefaultAssay(seu_obj) <- "Metab"
    if (!("samID" %in% colnames(seu_obj[[]]))) seu_obj$samID = seu_obj$patientID
    if (tumor == "PRAD") seu_obj$samID <- seu_obj$samID2
    
    epi <- subset(seu_obj, celltype == "Malignant" & TorN == "T")
    epi.list <- SplitObject(epi, split.by = "samID")
    nCells = sapply(epi.list, function(x) dim(x)[2])
    epi.list <- epi.list[nCells > 20]
    
    epi.list <- lapply(epi.list, function(seu) {
        DefaultAssay(seu) = "Metab"
        seu <- FindVariableFeatures(seu, verbose = FALSE)
    })
    varFeatures.list <- lapply(epi.list, function(seu) {
        res = seu@assays$Metab@meta.features
        idx = order(res$vst.variance.standardized, decreasing = TRUE)
        return(head(rownames(res)[idx], 10))
    })
    varFeatures <- data.frame(table(unlist(varFeatures.list)))
    varFeatures$Percent <- varFeatures$Freq/length(varFeatures.list)
    varFeatures$dataset <- tumor
    intrasam_feature[[tumor]] <- varFeatures
}
a = intrasam_feature
intrasam_feature <- do.call(rbind, intrasam_feature)
write.table(intrasam_feature, file = "res/3_epi/IntraSample/intrasam_feature.tsv", sep = "\t", row.names = F, quote = F)


intrasam_feature <- read.delim("res/3_epi/IntraSample/intrasam_feature.tsv")
intrasam_feature <- subset(intrasam_feature, Percent > 0.2)
intrasam_feature$Var1 <- gaude_trans[intrasam_feature$Var1, ]$updated

# plot 
label_data <- intrasam_feature %>%
    group_by(Var1) %>%
    summarise(n = n(), stats = sum(Percent))
label_data <- label_data[order(label_data$n, label_data$stats, decreasing = TRUE), ]
intrasam_feature$Var1 <- factor(intrasam_feature$Var1, levels = label_data$Var1)
label_data$id <- 1:nrow(label_data)
angle <- 90 - 360 * (label_data$id - 0.5) / nrow(label_data)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)


tmpfile <- paste0("plot/", dir.used, "PanCan_intra_sample_stacked_change_color.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 6, height = 6)
p <- ggplot(intrasam_feature, aes(x = Var1, y = Percent, fill = dataset)) + 
    geom_bar(stat = "identity", alpha = 1) + 
    ylim(-2, max(label_data$stats) + 18) + 
    theme_minimal() +
    scale_fill_manual(values = dataset_color) + 
    theme(
        legend.position = "right", 
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        text = tmp.text
    ) +
    coord_polar() + 
    geom_text(data = label_data, aes(x = id, y = stats + 1, label = Var1, hjust = hjust), 
              color = "black", size = 2.82, angle = label_data$angle, inherit.aes = FALSE)
p
dev.off()
