library(Seurat)
options(Seurat.object.assay.version = "v5")
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
options(future.globals.maxSize = 1e+09)

# 0. prepare file for PADC_xxxx ==========================================
# tmpobj <- Read10X("data/NMF_valid/ST/ST-PDAC_HT231P1_NG_2022_Zhou/")
# DropletUtils::write10xCounts("data/NMF_valid/ST/ST-PDAC_HT231P1_NG_2022_Zhou/filtered_feature_bc_matrix.h5",
#                              x = tmpobj)
# 
# tmpobj <- Read10X("data/NMF_valid/ST/ST-PDAC_HT259P1_NG_2022_Zhou/")
# DropletUtils::write10xCounts("data/NMF_valid/ST/ST-PDAC_HT259P1_NG_2022_Zhou/filtered_feature_bc_matrix.h5",
#                              x = tmpobj)
# 
# # I. run in server ===============================================
# MP_list <- readRDS("res/02_NMF/MPs/MP_list.RDS")$Malignant
# for (sam in dir("data/NMF_valid/ST/")) {
#     stobj <- Load10X_Spatial(paste0("data/NMF_valid/ST/", sam))
#     
#     # calculate MMP score
#     stobj <- NormalizeData(stobj, scale.factor = 1e5, verbose = F)
#     expr <- GetAssayData(stobj, assay = "Spatial", layer = "data") / log(2) ## transform to log2(CPM/10 + 1)
#     rownames(expr) <- scMetab::update_symbols(rownames(expr))$updated
#     MMP_score = scalop::sigScores(as.matrix(expr), sigs = MP_list, conserved.genes = 0.2)
#     stobj <- AddMetaData(stobj, MMP_score)
#     
#     # SpatialFeaturePlot(stobj, features = "MMP9")
#             
#     stobj <- stobj %>% 
#         SCTransform(assay = "Spatial", verbose = FALSE) %>%
#         RunPCA(assay = "SCT", verbose = FALSE) %>%
#         FindNeighbors(reduction = "pca", dims = 1:30) %>%
#         RunUMAP(reduction = "pca", dims = 1:30) %>% 
#         FindClusters(verbose = FALSE, resolution = c(2:8)/10)
#     saveRDS(object = stobj, file = paste0("res/02_NMF/MPs/ST/", sam, ".RDS"))
# }

# II. plot MMP scores ==============================================================================
MP_list <- readRDS("res/02_NMF/MPs/MP_list.RDS")$Malignant
sam = "ST-PDAC_A1"

for (sam in dir("data/NMF_valid/ST/", pattern = "ST-PDAC_[ABC]1")) {
    stobj <- readRDS(file = paste0("res/02_NMF/MPs/ST/", sam, ".RDS"))
    
    SpatialDimPlot(stobj, group.by = paste0("SCT_snn_res.", c(2:8)/10))
    
    
    anno <- read.csv(paste0("data/NMF_valid/ST/", sam, "/my_anno.csv"))
    #all(anno$Barcode == colnames(stobj))
    stobj$Group = anno$my.category[match(colnames(stobj), anno$Barcode)]
    stobj$Group = ifelse(stobj$Group == "", NA, 
                         ifelse(stobj$Group %in% c("HG PanIN", "LG PanIN"), "PanIN", 
                                ifelse(stobj$Group == "PDAC", "PAAD", stobj$Group)))
    group_color = celltype_color[1:3]
    names(group_color) = c("PAAD", "normal pancreas", "PanIN")
    
    p1 = SpatialDimPlot(stobj, group.by = "Group", pt.size.factor = 0) + NoLegend()
    p2 = SpatialDimPlot(stobj, group.by = "Group", cols = group_color)
    p3 = SpatialFeaturePlot(stobj, features = "MMP9")
    
    df = stobj[[c("MMP9", "Group")]] %>% filter(!is.na(Group))
    p4 = ggplot(df, aes(x = Group, MMP9, fill = Group)) +
        geom_boxplot() + 
        scale_fill_manual(values = group_color) + 
        theme_classic() + 
        stat_compare_means(comparisons = list(c("PAAD", "PanIN"), c("normal pancreas", "PanIN")))
    ggsave(p1+p2+p3+p4, file = paste0("plot/02_NMF/ST/", sam, "_combined.pdf"), width = 7, height = 7)
    
    p5 = SpatialFeaturePlot(stobj, features = colnames(MP_list))
    ggsave(p5, file = paste0("plot/02_NMF/ST/", sam, "_All_MMPs.pdf"), width = 16, height = 16)
}



for (sam in dir("data/NMF_valid/ST/", pattern = "ST-PDAC_HT", )) {
    stobj <- readRDS(file = paste0("res/02_NMF/MPs/ST/", sam, ".RDS"))
    
    SpatialDimPlot(stobj, group.by = paste0("SCT_snn_res.", c(2:8)/10))
    
    p1 = SpatialDimPlot(stobj, pt.size.factor = 0) + NoLegend()
    p2 = SpatialDimPlot(stobj, group.by = "SCT_snn_res.0.6")
    p3 = SpatialFeaturePlot(stobj, features = "MMP9")
    ggsave(p1+p2+p3, file = paste0("plot/02_NMF/ST/", sam, "_combined.pdf"), width = 10.5, height = 3.5)
}


for (sam in dir("data/NMF_valid/ST/", pattern = "ST-colon", )) {
    stobj <- readRDS(file = paste0("res/02_NMF/MPs/ST/", sam, ".RDS"))
    
    # SpatialDimPlot(stobj, group.by = paste0("SCT_snn_res.", c(2:8)/10))
    
    p1 = SpatialDimPlot(stobj, pt.size.factor = 0) + NoLegend()
    p3 = SpatialFeaturePlot(stobj, features = "MMP7")
    ggsave(p1+p3, file = paste0("plot/02_NMF/ST/", sam, "_combined.pdf"), width = 7, height = 3.5)
    
    p5 = SpatialDimPlot(stobj, group.by = paste0("SCT_snn_res.", c(2:8)/10))
    ggsave(p5, file = paste0("plot/02_NMF/ST/", sam, "_All_MMPs.pdf"), width = 12, height = 12)
}







# anno <- read.csv("/mnt/e/G_Bak/DD/STdata/data_for_analysis/spots_annotation/1622021_case1_annotation.csv")
# stobj$Group <- anno$my.category

SpatialFeaturePlot(stobj, features = "nCount_Spatial")
stobj <- SCTransform(stobj, assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(stobj, features = "TFF2")


stobj <- stobj %>% 
    RunPCA(assay = "SCT", verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    RunUMAP(reduction = "pca", dims = 1:30)
stobj <- FindClusters(stobj, verbose = FALSE, resolution = 0.3) 



# expr


stobj <- AddModuleScore(stobj, MP_list, assay = "SCT")


stobj <- AddMetaData(stobj, res)
SpatialFeaturePlot(stobj, features = paste0("Cluster", 1:15)) 
SpatialDimPlot(stobj,group.by = "Group")

p1 <- DimPlot(stobj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(stobj, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(stobj, label = TRUE, label.size = 3, group.by = "Group")
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
de_markers <- FindAllMarkers(stobj, only.pos = T)