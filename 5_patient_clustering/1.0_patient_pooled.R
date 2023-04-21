# ==================================================================================================
# Author: Zhou Zhe
# Program: prepare data for SCENIC
#          Datasets: for all datasets
#          Samples: ..
#          TorN: only Tumor samples
#          Celltype: All celltype included
# Version: 1.0
# Date: Jul 7, 2022
# ==================================================================================================

# 0. INIT -----------------------------------------------------------------------------------------
rm(list=ls())


library(tidyverse)
library(Seurat)
library(SCopeLoomR)
library(scMetab)
library(ComplexHeatmap)
library(GSVA)
library(parallel)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(dendextend)
library(circlize)

dir.used <- "res/6_patient_clustering/"
file.prefix <- "patient_celltype"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
seed.use = 2021

# I. Pool cells -----------------------------------------------------------------------------------
## 1. Pool cells by sample x celltype ======
pooled_list <- lapply(1:nrow(dataset), function(i) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    # dir.create(dir.used, recursive = T)
    # sink(paste0(dir.used, "log/", file.prefix, tumor, ".log"))
    
    seu_obj <- readRDS(paste0("data/rds/PAS", tumor, ".TpN.PAS.unintegrated.rds"))
    seu_obj <- subset(seu_obj, TorN == "T")
    # seu_obj <- readRDS(paste0("../", dataset$Directory[i], "/data/rds/", dataset$Prefix[i], ".merged.rds"))
    
    if (!("samID" %in% colnames(seu_obj[[]]))) seu_obj$samID = seu_obj$patientID
    if (tumor == "PRAD") seu_obj$samID <- seu_obj$samID2
    
    # 倾向于使用samID进行pool
    
    seu_obj$samID_celltype = paste(tumor, seu_obj$samID, seu_obj$celltype, sep = "_")
    meta_info <- seu_obj[[c("samID_celltype", "samID", "celltype", "patientID")]]
    meta_info <- meta_info %>%
        group_by(samID_celltype) %>%
        summarise(n = n(), celltype = unique(celltype), samID = unique(samID), patientID = unique(patientID)) %>%
        as.data.frame()
    meta_info$tumor = tumor
    meta_info <- meta_info[order(meta_info$celltype, meta_info$n, decreasing = T), ]
    rownames(meta_info) <- meta_info$samID_celltype
    # meta_info <- subset(meta_info, n >= 20)
    
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj <- NormalizeData(seu_obj)
    pooled_obj <- AverageExpression(seu_obj, assays = c("RNA", "Metab"), slot = c("data", "counts"), 
                                    group.by = "samID_celltype", return.seurat = T)
    pooled_obj <- pooled_obj[, colnames(pooled_obj) %in% rownames(meta_info)]
    pooled_obj <- AddMetaData(pooled_obj, meta_info)
    return(pooled_obj)
})

pooled.merged <- merge(pooled_list[[1]], pooled_list[-1])
pooled.merged@assays$Metab@data <- GetAssayData(pooled.merged, assay = "Metab", slot = "counts")
fastSave::saveRDS.pigz(pooled.merged, file = "data/rds/pseudo_sam_celltype_tmp.rds", n.cores = 10)

# clean and save data
DefaultAssay(pooled.merged) <- "Metab"
pooled.merged <- ScaleData(pooled.merged)

DefaultAssay(pooled.merged) <- "RNA"
metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(pooled.merged))
pooled.merged <- NormalizeData(pooled.merged)
VariableFeatures(pooled.merged) <- metab.gene
pooled.merged <- ScaleData(pooled.merged)

rna_counts <- GetAssayData(pooled.merged, assay = "RNA", slot = "counts")[metab.gene, ]
pooled.merged$nFeature_metab <- colSums(rna_counts > 0)
pooled.merged$nCount_metab <- colSums(rna_counts)
fastSave::saveRDS.pigz(pooled.merged, file = "data/rds/pseudo_sam_celltype.rds", n.cores = 10)


## 2. pooled cells by sample ======
pooled_list_sam <- lapply(1:nrow(dataset), function(i) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    # dir.create(dir.used, recursive = T)
    # sink(paste0(dir.used, "log/", file.prefix, tumor, ".log"))
    
    seu_obj <- readRDS(paste0("data/rds/", tumor, ".TpN.PAS.unintegrated.rds"))
    seu_obj <- subset(seu_obj, TorN == "T")
    
    if (!("samID" %in% colnames(seu_obj[[]]))) seu_obj$samID = seu_obj$patientID
    if (tumor == "PRAD") seu_obj$samID <- seu_obj$samID2
    
    # 倾向于使用samID进行pool
    
    # seu_obj$samID_celltype = paste(tumor, seu_obj$samID, seu_obj$celltype, sep = "_")
    meta_info <- seu_obj[[c("samID", "patientID")]]
    meta_info <- meta_info %>%
        group_by(samID) %>%
        summarise(n = n(), patientID = unique(patientID)) %>%
        as.data.frame()
    meta_info$dataset = tumor
    meta_info$cancerType = sub("_.+", "", tumor)
    meta_info <- meta_info[order(meta_info$n, decreasing = T), ]
    rownames(meta_info) <- meta_info$samID
    # meta_info <- subset(meta_info, n >= 20)
    
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj <- NormalizeData(seu_obj)
    pooled_obj <- AverageExpression(seu_obj, assays = c("RNA", "Metab"), slot = c("data", "counts"), 
                                    group.by = "samID", return.seurat = T)
    pooled_obj <- pooled_obj[, colnames(pooled_obj) %in% rownames(meta_info)]
    pooled_obj <- AddMetaData(pooled_obj, meta_info)
    return(pooled_obj)
})
pooled.merged.sam <- merge(pooled_list_sam[[1]], pooled_list_sam[-1])
pooled.merged.sam@assays$Metab@data <- GetAssayData(pooled.merged.sam, assay = "Metab", slot = "counts")

# clean and save data
DefaultAssay(pooled.merged.sam) <- "RNA"
pooled.merged.sam <- NormalizeData(pooled.merged.sam)
fastSave::saveRDS.pigz(pooled.merged.sam, file = "data/rds/pseudo_sam.rds", n.cores = 10)



# II. GSVA ----------------------------------------------------------------------------------------
metab.gene.list <- split(gaudeMetabDf$GeneSymbol, f = gaudeMetabDf$Pathway)
metab.gene.list <- metab.gene.list[sapply(metab.gene.list, length) >=5]

rna_data <- GetAssayData(pooled.merged, assay = "RNA", slot = "data")
gene_info <- data.frame(nCells = rowSums(rna_data > 0), nCounts = rowSums(rna_data))
rna_data <- rna_data[gene_info$nCells > 0, ]

metab.gsva <- gsva(as.matrix(rna_data), metab.gene.list, method = "gsva", kcdf = "Gaussian", mx.diff = T, parallel.sz = 10)
fastSave::saveRDS.pigz(metab.gsva, file = "data/rds/pseudo_sam_celltype_metab_gsva.rds")



# III. Cluster Pathway -----------------------------------------------------------------------------
pooled.merged <- readRDS("data/rds/pseudo_sam_celltype.rds")
pooled.merged$cancerType = pooled.merged$orig.ident
pooled.merged$dataset = pooled.merged$tumor
metab.gsva <- readRDS("data/rds/pseudo_sam_celltype_metab_gsva.rds")
pooled.merged <- subset(pooled.merged, n >= 20)
metab.gsva <- metab.gsva[, colnames(pooled.merged)]

#### consensusClusterPlus
title = "plot/6_patient_clustering/Pathway_Consensus"
consensus_path <- ConsensusClusterPlus(t(metab.gsva), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, seed = seed.use, 
                                     title = title, clusterAlg = "pam", distance = "spearman", plot = "pdf")
icl_path <- calcICL(consensus_path, title = title, plot = "pdf")
save(consensus_path, icl_path, file = "plot/6_patient_clustering/Pathway_Consensus/consensus_res.rda")

#### according to the elbow plot, choose 5 as the number of pathway modules
k = 5
pathway_dend <- consensus_path[[k]]$consensusTree
pathway_class <- data.frame(Module = paste0("Module_", consensus_path[[k]]$consensusClass), 
                            row.names = names(consensus_path[[k]]$consensusClass))
pathway_dend$labels <- gaude_trans[rownames(pathway_class), ]$updated

#### circlize dendrogram plot
hc <- as.dendrogram(pathway_dend) %>%
    color_branches(k = 5) %>%
    color_labels(k = 5)
pdf(file.path(title, "pathway_circlize_dend_resize.pdf"), width = 6, height = 6)
par(ps = 6)
circlize_dendrogram(hc, dend_track_height = 0.3, labels_track_height = 0.63)
dev.off()

#### use color coincide with dendrogram for Module
dend_color <- colorspace::rainbow_hcl(5, c=90, l=50)
names(dend_color) <- unique(pathway_class[pathway_dend$order, ])

# path_cor <- cor(t(metab.gsva), method = "pearson")
pdf(file.path(title, "pathway_consensus.pdf"), width = 6, height = 4)
tmp.text <- gpar(fontfamily="sans", fontsize = 16)
Heatmap(consensus_path[[4]]$consensusMatrix, col = circlize::colorRamp2(c(0, 1), c("white", brewer.pal(8, "Accent")[5])), 
        top_annotation = columnAnnotation(df = pathway_class, col = list(Module = dend_color), 
                                          simple_anno_size = unit(2, "mm"), show_annotation_name = F, 
                                          annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text)), 
        # left_annotation = rowAnnotation(df = pathway_class, col = list(Module = dend_color), 
        #                                 simple_anno_size = unit(2, "mm"), show_legend = F, show_annotation_name = F, 
        #                                 annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text)), 
        show_column_dend = T, show_row_dend = F, 
        name = "consensus", 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        cluster_rows = pathway_dend, cluster_columns = pathway_dend,
        row_names_gp = tmp.text, column_names_gp = tmp.text)
dev.off()


# IV. Cluster Sample ------------------------------------------------------------------------------
# heatmap function
metab.heatmap <- function(plot_data, 
                          pdffile = NULL, 
                          show_row_names = T, 
                          column_km = 1, 
                          column_split = NULL) {
    sam_Info <- pooled.merged[[c("celltype", "cancerType", "dataset")]]
    cancerType_color = dataset_color[c(1, 4, 5, 8:10)]
    
    sams_used <- colnames(plot_data)
    sam_Info <- sam_Info[sams_used, ]
    
    
    pdf(pdffile, width = 20, height = 12)
    col_ha = columnAnnotation(df = sam_Info, col = list(celltype = celltype_color,
                                                        cancerType = cancerType_color, 
                                                        dataset = dataset_color))
    withr::local_seed(2021)
    p <- Heatmap(plot_data, top_annotation = col_ha,
                 colorRamp2(c(-2, 0, 2), brewer.pal(11, "RdBu")[c(10, 6, 2)]),
                 show_column_names = F, show_row_names = show_row_names, 
                 right_annotation = rowAnnotation(df = pathway_class, col = list(Module = dend_color), 
                                                 simple_anno_size = unit(2, "mm"), show_legend = F, show_annotation_name = F, 
                                                 annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text)), 
                 clustering_distance_columns = "spearman", 
                 cluster_rows = pathway_dend, cluster_columns = T, 
                 row_split = 5, column_split = column_split, 
                 column_km = column_km, column_km_repeats = 100)
    p <- draw(p)
    
    dev.off()
    
    return(p)
}



samClu_path = "plot/6_patient_clustering/samClu/"
plot_data <- t(scale(t(metab.gsva)))
metab.heatmap(plot_data = plot_data, column_km = 8,
              pdffile = file.path(samClu_path, "Metab_gsva_clustering_ALL_change_color.pdf"))


for (used_celltype in c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")) {
    sams_used <- colnames(pooled.merged)[pooled.merged$celltype == used_celltype]
    plot_res = metab.heatmap(plot_data = plot_data[, sams_used], column_km = 2, 
                  pdffile = paste0(samClu_path, "Metab_gsva_clustering_", used_celltype, "_change_color_ymp.pdf"))
    
    # sams <- plot_res@ht_list[[1]]@column_names_param$labels
    # sam_clu <- reshape2::melt(column_order(plot_res))
    # sam_clu$name = sams[sam_clu$value]
    # sam_clu <- cbind(sam_clu, pooled.merged[[c("celltype", "cancerType", "dataset", "samID")]][sam_clu$name, ])
    # sam_clu$cluster <- paste0(sub(" ", "-", sam_clu$celltype), "_C", sam_clu$L1)
    # write.csv(sam_clu[c(-1, -2)], file = paste0("res/6_patient_clustering/samClu/", used_celltype, "_Cluster.csv"), 
    #           quote = FALSE, row.names = F)
    # 
    # print(table(sam_clu$cluster, sam_clu$dataset))
}


# V. Cluster DE -----------------------------------------------------------------------------------
pooled.merged <- readRDS("data/rds/pseudo_sam_celltype.rds")
pooled.merged.sam <- readRDS("data/rds/pseudo_sam.rds")

for (used_celltype in c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")) {
    sam_clu <- read.csv(paste0("res/6_patient_clustering/", used_celltype, "_Cluster.csv"))
    
    Idents(pooled.merged) <- sam_clu$cluster[match(colnames(pooled.merged), sam_clu$name)]
    table(Idents(pooled.merged), useNA = "ifany")
    markers_celltype <- FindMarkers(pooled.merged, slot = "data", 
                                    ident.1 = paste0(sub(" ", "-", used_celltype), "_C1"),
                                    ident.2 = paste0(sub(" ", "-", used_celltype), "_C2"))
    table(markers_celltype$avg_log2FC > 0)
    write.csv(markers_celltype, file = paste0("res/6_patient_clustering/Cluster_DEG/", used_celltype, "_celltype_DEG.csv"))
    
    Idents(pooled.merged.sam) <- sam_clu$cluster[match(colnames(pooled.merged.sam), sam_clu$samID)]
    table(Idents(pooled.merged.sam), useNA = "ifany")
    markers_sam <- FindMarkers(pooled.merged.sam, slot = "data", 
                               ident.1 = paste0(sub(" ", "-", used_celltype), "_C1"),
                               ident.2 = paste0(sub(" ", "-", used_celltype), "_C2"))
    table(markers_sam$avg_log2FC > 0)
    write.csv(markers_sam, file = paste0("res/6_patient_clustering/Cluster_DEG/", used_celltype, "_sample_DEG.csv"))    
}


# VI. Cluster DE Pathway Analysis -----------------------------------------------------------------
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
celltypes <- c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")
cluster_sigs <- lapply(celltypes, function(used_celltype) {
    markers_celltype <- read.csv(file = paste0("res/6_patient_clustering/Cluster_DEG/", used_celltype, "_celltype_DEG.csv"))
    markers_celltype <- subset(markers_celltype, p_val_adj < 0.01)
    markers_celltype <- markers_celltype %>%
        group_by(avg_log2FC > 0) %>%
        slice_max(order_by = abs(avg_log2FC), n = 50)
    
    markers_sam <- read.csv(file = paste0("res/6_patient_clustering/Cluster_DEG/", used_celltype, "_sample_DEG.csv"))
    markers_sam <- subset(markers_sam, p_val_adj < 0.05)
    markers_sam <- markers_sam %>%
        group_by(avg_log2FC > 0) %>%
        slice_max(order_by = abs(avg_log2FC), n = 50)    
    
    res = list(celltype_C1 = markers_celltype$X[markers_celltype$avg_log2FC > 0], 
               celltype_C2 = markers_celltype$X[markers_celltype$avg_log2FC < 0], 
               sam_C1 = markers_sam$X[markers_sam$avg_log2FC > 0], 
               sam_C2 = markers_sam$X[markers_sam$avg_log2FC < 0])
    names(res) = paste0(sub(" ", "-", used_celltype), "_", names(res))
    return(res)
})
cluster_sigs <- do.call(c, cluster_sigs)

library(clusterProfiler)
library(org.Hs.eg.db)
gl <- cluster_sigs$Myeloid_sam_C2
egoBP <- enrichGO(gene = gl, 
                  OrgDb = org.Hs.eg.db, 
                  keyType = "SYMBOL", 
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1)









