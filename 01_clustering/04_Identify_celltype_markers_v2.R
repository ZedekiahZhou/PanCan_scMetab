# ====================================================================
# Author: Zhou Zhe
# Function: Plot celltype specific metab for all cancers
# Version: 2.0
# History:
#   1.0: Only Pathway 
#   2.0: Add TF
# Date: Jul 18, 2022
# ====================================================================

rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyverse)
library(scMetab)
library(patchwork)
library(Seurat)
library(GSEABase)
library(grid)

dir.used <- "1_clustering/PAS_markers/"
file.prefix <- "Celltype_"

pseudo <- 1e-50
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

gaude.list <- split(gaudeMetabDf$GeneSymbol, f = gaudeMetabDf$Pathway)

# I. Prepare data ===========================================================================================
## 1. Gene markers ------
tmpfile <- paste0("plot/", dir.used, file.prefix, "Gene_split.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 6, height = 6)
gene_markers.list <- lapply(rownames(dataset), function(tumor) {
    gene_markers <- read.delim(paste0("res/1_clustering/PAS_clustering/", tumor, "_Celltype_Gene_Markers.tsv"))
    gene_markers$log10padj <- (-1) * sign(gene_markers$avg_log2FC) * log10(gene_markers$p_val_adj)
    gene_markers$tumor <- tumor
    gene_markers$log10padj_pseudo <- ifelse(gene_markers$log10padj > 50, 50, 
                                            ifelse(gene_markers$log10padj < -50, -50, gene_markers$log10padj))
    gene_markers <- subset(gene_markers, (p_val_adj < 0.01) & (avg_log2FC > 0.1) & ((pct.1 > 0.1) | (pct.2 > 0.1)))
    gene_markers <- gene_markers[order(gene_markers$cluster, -gene_markers$log10padj), ]
    
    
    top10 <- gene_markers %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = log10padj)
    top10 <- top10
    top10 <- top10[order(top10$cluster, top10$p_val_adj), ]
    top10$gene <- factor(top10$gene, levels = unique(top10$gene))
    mycolor <- brewer.pal(11, "RdBu")[c(10, 6, 2)]
    print(ggplot(top10, aes(x = cluster, y = gene, size = abs(avg_log2FC), color = log10padj_pseudo)) + 
              geom_point() + 
              # scale_color_distiller(palette = "RdBu") + 
              ggtitle(tumor) +
              scale_size(limits = c(0, max(abs(top10$avg_log2FC))), range = c(0.5, 2.5)) +
              scale_color_gradient2(low = mycolor[1], high = mycolor[3], mid = mycolor[2]) + 
              theme_grey() + 
              theme(text=tmp.text, axis.text.x = element_text(angle = 45, hjust = 1), 
                    axis.line=element_line(size = .3, colour="black")))
    # gene_markers.filtered <- subset(gene_markers, (p_val_adj < 0.01) & (avg_log2FC > 0.25) & ((pct.1 > 0.1) | (pct.2 > 0.1)))
    return(gene_markers)
})
dev.off()
names(gene_markers.list) <- rownames(dataset)
celltype_gene_markers <- do.call(rbind, gene_markers.list)
write.table(celltype_gene_markers, file = paste0("res/", dir.used, file.prefix, "Gene_Markers_All_Datasets.tsv"), 
            quote = F, sep = "\t", row.names = F)

## 2. Metabolic markers ------
tmpfile <- paste0("plot/", dir.used, file.prefix, "metab_split.pdf")
tmp.text <- element_text(family="sans", size=21)
pdf(tmpfile, width = 12, height = 13)
metab_markers.list <- lapply(rownames(dataset), function(tumor) {
    metab_markers <- read.delim(paste0("res/1_clustering/PAS_clustering/", tumor, "_Celltype_Metab_Markers.tsv"))
    metab_markers$log10padj <- (-1) * sign(metab_markers$avg_log2FC) * log10(metab_markers$p_val_adj)
    metab_markers$log10padj_pseudo <- ifelse(metab_markers$log10padj > 50, 50, 
                                             ifelse(metab_markers$log10padj < -50, -50, metab_markers$log10padj))
    metab_markers$tumor <- tumor
    metab_markers <- subset(metab_markers, (p_val_adj < 0.01) & (abs(avg_log2FC) > 0.01) & ((pct.1 > 0.1) | (pct.2 > 0.1)))
    metab_markers <- metab_markers[order(metab_markers$cluster, -metab_markers$log10padj), ]
    metab_markers$gene <- factor(metab_markers$gene, levels = unique(metab_markers$gene))
    ### core enrichment
    metab_markers$core <- apply(metab_markers[, c("cluster", "gene")], 1, function(x) {
        genes = intersect(gaude.list[[x[2]]], gene_markers.list[[tumor]]$gene[gene_markers.list[[tumor]]$cluster == x[1]])
        return(paste0(genes, collapse = "/"))
    })
    
    mycolor <- brewer.pal(11, "RdBu")[c(10, 6, 2)]
    print(ggplot(metab_markers, aes(x = cluster, y = gene, size = abs(avg_log2FC), color = log10padj_pseudo)) + 
              geom_point() + 
              # scale_color_distiller(palette = "RdBu") + 
              ggtitle(tumor) +
              scale_size(limits = c(0, max(abs(metab_markers$avg_log2FC))), range = c(0.5, 4)) +
              scale_color_gradient2(low = mycolor[1], high = mycolor[3], mid = mycolor[2]) + 
              theme(text=tmp.text, axis.text.x = element_text(angle = 45, hjust = 1), 
                    axis.line=element_line(size = .3, colour="black")))
    return(metab_markers)
    
})
dev.off()
names(metab_markers.list) <- rownames(dataset)
celltype_metab_markers <- do.call(rbind, metab_markers.list)
write.table(celltype_metab_markers, file = paste0("res/", dir.used, file.prefix, "Metab_Markers_All_Datasets.tsv"), 
            quote = F, sep = "\t", row.names = F)

## 3. TF markers ------
tmpfile <- paste0("plot/", dir.used, file.prefix, "TF_split.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 6, height = 6)
TF_markers.list <- lapply(rownames(dataset), function(tumor) {
    TF_markers <- read.delim(paste0("res/1_clustering/PAS_clustering/", tumor, "_Celltype_TF_Markers.tsv"))
    TF_markers$log10padj <- (-1) * sign(TF_markers$avg_log2FC) * log10(TF_markers$p_val_adj)
    TF_markers$tumor <- tumor
    TF_markers$log10padj_pseudo <- ifelse(TF_markers$log10padj > 50, 50, 
                                            ifelse(TF_markers$log10padj < -50, -50, TF_markers$log10padj))
    TF_markers <- subset(TF_markers, (p_val_adj < 0.01) & (abs(avg_log2FC) > 0.01) & ((pct.1 > 0.1) | (pct.2 > 0.1)))
    TF_markers <- TF_markers[order(TF_markers$cluster, -TF_markers$log10padj), ]
    
    ### core enrichment
    regulons <- geneIds(getGmt(paste0("res/5_SCENIC/regulon/", tumor, "_reg_metab.gmt")))
    TF_markers$core <- apply(TF_markers[, c("cluster", "gene")], 1, function(x) {
        genes = intersect(regulons[[x[2]]], gene_markers.list[[tumor]]$gene[gene_markers.list[[tumor]]$cluster == x[1]])
        return(paste0(genes, collapse = "/"))
    })
    
    top10 <- TF_markers %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = log10padj)
    top10 <- top10
    top10 <- top10[order(top10$cluster, top10$p_val_adj), ]
    top10$gene <- factor(top10$gene, levels = unique(top10$gene))
    mycolor <- brewer.pal(11, "RdBu")[c(10, 6, 2)]
    print(ggplot(top10, aes(x = cluster, y = gene, size = abs(avg_log2FC), color = log10padj_pseudo)) + 
              geom_point() + 
              # scale_color_distiller(palette = "RdBu") + 
              ggtitle(tumor) +
              scale_size(limits = c(0, max(abs(top10$avg_log2FC))), range = c(0.5, 2.5)) +
              scale_color_gradient2(low = mycolor[1], high = mycolor[3], mid = mycolor[2]) + 
              theme_grey() + 
              theme(text=tmp.text, axis.text.x = element_text(angle = 45, hjust = 1), 
                    axis.line=element_line(size = .3, colour="black")))
    return(TF_markers)
})
dev.off()
names(TF_markers.list) <- rownames(dataset)
celltype_TF_markers <- do.call(rbind, TF_markers.list)
write.table(celltype_TF_markers, file = paste0("res/", dir.used, file.prefix, "TF_Markers_All_Datasets.tsv"), 
            quote = F, sep = "\t", row.names = F)



# II. Plot Pooled multiple cancer ======================================================================
## 1. Metabolic markers ---------
celltype_metab_markers <- read.delim(paste0("res/", dir.used, file.prefix, "Metab_Markers_All_Datasets.tsv"))
pooled_metab <- celltype_metab_markers
pooled_metab <- read.delim(paste0("res/", dir.used, file.prefix, "Metab_Markers_All_Datasets.tsv"))
pooled_metab <- pooled_metab[pooled_metab$avg_log2FC > 0, ]
used_celltypes <- c("Malignant", "Myeloid", "Fibroblasts", "Endothelial", "T cells", "B cells")
pooled_metab <- pooled_metab[pooled_metab$cluster %in% used_celltypes, ]

### Heatmap of number of datasets
mat <- as.data.frame.matrix(table(pooled_metab$gene, pooled_metab$cluster))
mat <- mat[, used_celltypes]
mat$max_celltype <- factor(apply(mat[1:6], 1, function(x) colnames(mat[1:6])[which.max(x)]), 
                           levels = used_celltypes)
mat$max <- apply(mat[1:6], 1, max)
mat <- mat[mat$max>0, ]
mat <- mat[order(mat$max_celltype, -mat$max), ]
rownames(mat) <- gaude_trans$updated[match(rownames(mat), gaude_trans$original)]

tmpfile <- paste0("plot/", dir.used, file.prefix, "MetabPath_pooled_All_10_datasets.pdf")
tmp.text <- gpar(fontfamily="sans", fontsize = 8)
pdf(tmpfile, width = 5, height = 9)
Heatmap(mat[1:6], circlize::colorRamp2(seq(0, 10, 2.5), brewer.pal(11, "RdBu")[6:2]), 
        rect_gp = gpar(col = "white"), 
        row_names_max_width = unit(6, "cm"),
        row_names_gp = tmp.text, column_names_gp = tmp.text, 
        row_title_gp = tmp.text, column_title_gp = tmp.text, 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        name = "Number of \n Cancers", 
        cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


## 2. Gene Markers list ----------
used_celltypes <- c("Malignant", "Myeloid", "Fibroblasts", "Endothelial", "T cells", "B cells")
celltype_gene_markers <- read.delim(paste0("res/", dir.used, file.prefix, "Gene_Markers_All_Datasets.tsv"))

### pie plot of conserved cell type-specific metabolic genes
n_all <- length(unique(gaudeMetabDf$GeneSymbol))
n_cts <- length(unique(celltype_gene_markers$gene))
prob_cts <- n_cts/n_all

df_con <- celltype_gene_markers %>% 
    group_by(cluster, gene) %>%
    summarise(n_con = n())
n_con <- length(unique(df_con$gene[df_con$n_con >= 8]))

plot_data <- data.frame(group = c("Non Sig", "Sig - Con", "Con"), 
                        value = c(n_all-n_cts, n_cts-n_con, n_con))
ggpubr::ggpie(plot_data, "value", label = "group")
ggsave(paste0("plot/", dir.used, "Conserve_metab_genes_pie.pdf"))

### Heatmap of number of datasets
pooled_gene <- subset(celltype_gene_markers, cluster %in% used_celltypes & avg_log2FC > 0.25)
mat <- as.data.frame.matrix(table(pooled_gene$gene, pooled_gene$cluster))
used_gene <- sapply(mat, function(x) rownames(mat)[order(x, decreasing = T)][1:5])
used_gene <- unique(as.vector(used_gene))

mat <- mat[used_gene, used_celltypes]
mat$max_celltype <- factor(apply(mat[1:6], 1, function(x) colnames(mat[1:6])[which.max(x)]), 
                           levels = used_celltypes)
mat$max <- apply(mat[1:6], 1, max)
mat <- mat[order(mat$max_celltype, -mat$max), ]

tmpfile <- paste0("plot/", dir.used, file.prefix, "MetabGene_pooled_All_10_datasets.pdf")
tmp.text <- gpar(fontfamily="sans", fontsize = 8)
pdf(tmpfile, width = 6.2, height = 1.2)
Heatmap(t(mat[1:6]), circlize::colorRamp2(seq(0, 10, 2.5), brewer.pal(11, "RdBu")[6:2]), 
        rect_gp = gpar(col = "white"), 
        row_names_max_width = unit(6, "cm"),
        row_names_gp = tmp.text, column_names_gp = tmp.text, 
        column_names_rot = 45, 
        row_title_gp = tmp.text, column_title_gp = tmp.text, 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        name = "Number of \n Cancers", 
        cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

### Heatmap of average log2FC in multiple dataset
pooled_gene <- subset(celltype_gene_markers, cluster %in% used_celltypes)
gene_summary <- lapply(split(pooled_gene, f = pooled_gene$cluster), function(x) {
    x <- x %>% 
        group_by(gene) %>% 
        summarise(mean_avg_log2FC = sum(avg_log2FC)/10, 
                  median_log2FC = median(avg_log2FC), 
                  max_log2FC = max(avg_log2FC), 
                  n = n()) %>%
        arrange(desc(n), desc(mean_avg_log2FC))
    return(x)
})
top_genes <- unique(unlist(lapply(gene_summary[used_celltypes], function(x) x$gene[1:5])))
write.csv(top_genes, "res/1_clustering/PAS_markers/Celltype_top_genes.csv", 
          row.names = F, col.names = F, quote = F)
# top_genes <- unique(unlist(top_genes[used_celltypes]))
pooled_gene <- pooled_gene[pooled_gene$gene %in% top_genes, ] 
pooled_gene$cluster <- factor(pooled_gene$cluster, levels = used_celltypes)
pooled_gene$gene <- factor(pooled_gene$gene, levels = top_genes)
mat <- reshape2::dcast(pooled_gene, gene ~ cluster + tumor, value.var = "avg_log2FC")
mat <- column_to_rownames(mat, var = "gene")
mat[is.na(mat)] <- 0

colInfo <- data.frame(Celltype = sub("_.+$", "", colnames(mat)), 
                      Dataset = sub("^[^_]+_", "", colnames(mat)), 
                      row.names = colnames(mat))
tmpfile <- paste0("plot/", dir.used, file.prefix, "MetabGene_Per_10_datasets_change_color.pdf")
tmp.text <- gpar(fontfamily="sans", fontsize = 8)
pdf(tmpfile, width = 9, height = 3)
colAnno <- columnAnnotation(df = colInfo, col = list(Celltype = celltype_color, Dataset = dataset_color), 
                            gp = gpar(col = "white"), annotation_label = c("Cell Type", "Dataset"), 
                            annotation_name_gp = tmp.text, 
                            simple_anno_size = unit(0.2, "cm"), 
                            annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text))
Heatmap(mat, circlize::colorRamp2(c(0, 3), brewer.pal(11, "RdBu")[c(6, 2)]), 
        top_annotation = colAnno, rect_gp = gpar(col = "white"), 
        cluster_columns = F, cluster_rows = F, 
        show_column_names = F, row_names_gp = tmp.text, 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        name = "Average log2FC")
dev.off()


#### plot genes counts
pooled_gene <- subset(celltype_gene_markers, cluster %in% used_celltypes & avg_log2FC > 0.25)
tmpfile <- paste0("plot/", dir.used, file.prefix, "Gene_counts.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 6, height = 2)
ggplot(pooled_gene, mapping = aes(x = tumor, fill = cluster)) +
    geom_bar(position = "dodge") + 
    scale_fill_manual(values = celltype_color[names(celltype_color) %in% used_celltypes]) + 
    theme(text = tmp.text, axis.text = tmp.text, legend.text = tmp.text, 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.line=element_line(size = .3, colour="black"),
          legend.key.size = unit(3, "mm"))
dev.off()


## 3. TF Markers --------
### Heatmap of number of datasets
pooled_TF <- subset(celltype_TF_markers, cluster %in% used_celltypes & avg_log2FC > 0)
mat <- as.data.frame.matrix(table(pooled_TF$gene, pooled_TF$cluster))
used_TF <- sapply(mat, function(x) rownames(mat)[order(x, decreasing = T)][1:5])
used_TF <- unique(as.vector(used_TF))

mat <- mat[used_TF, used_celltypes]
mat$max_celltype <- factor(apply(mat[1:6], 1, function(x) colnames(mat[1:6])[which.max(x)]), 
                           levels = used_celltypes)
mat$max <- apply(mat[1:6], 1, max)
mat <- mat[order(mat$max_celltype, -mat$max), ]

tmpfile <- paste0("plot/", dir.used, file.prefix, "TF_pooled_All_10_datasets.pdf")
tmp.text <- gpar(fontfamily="sans", fontsize = 8)
pdf(tmpfile, width = 6.2, height = 1.2)
Heatmap(t(mat[1:6]), circlize::colorRamp2(seq(0, 10, 2.5), brewer.pal(11, "RdBu")[6:2]), 
        rect_gp = gpar(col = "white"), 
        row_names_max_width = unit(6, "cm"),
        row_names_gp = tmp.text, column_names_gp = tmp.text, 
        column_names_rot = 45, 
        row_title_gp = tmp.text, column_title_gp = tmp.text, 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        name = "Number of \n Cancers", 
        cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


# III. VlnPlot ====================================================================================
## 1. Metabolic Pathway ---------
seu.list <- lapply(c(1, 4, 7:10), function(i) {
    tumor <- dataset$DataSets[i]
    print(tumor)
    seu_obj <- readRDS(paste0("data/rds/PAS/", tumor, ".TpN.PAS.unintegrated.rds"))
    seu_obj[, seu_obj$TorN == "T"]
    gc()
    return(seu_obj)
})
names(seu.list) <- dataset$DataSets[c(1, 4, 7:10)]

# Plot Vlnplot
used_path <- c("Glycolysis and Gluconeogenesis", "ROS Detoxification", "Glutathione Metabolism", 
               "Pentose Phosphate Pathway", "Pterin Biosynthesis", "Eicosanoid Metabolism", 
              "Oxidative Phosphorylation", "Porphyrin and Heme Metabolism", "Transport, Lysosomal", 
              "Sphingolipid Metabolism", "Glutamate metabolism", "NAD Metabolism", 
              "Phosphoinositide Signalling", "Lysine Metabolism")
tmpfile <- paste0("plot/", dir.used, file.prefix, "PAS_VlnPlot.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 12, height = 9)
tmp <- lapply(names(seu.list), function(tumor) {
    seu <- seu.list[[tumor]]
    DefaultAssay(seu) <- "Metab"
    features = sapply(used_path, grep, x = rownames(seu), fixed = T, value = T)
    print(Seurat::VlnPlot(seu, features = features,
            group.by = "celltype", cols = celltype_color, pt.size = 0) &
        geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
        plot_annotation(title = tumor) &
        theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text))
})
dev.off()


used_path <- c("Glycolysis and Gluconeogenesis", "Oxidative Phosphorylation", "Glutamate metabolism")
used_celltypes <- c("Malignant", "Myeloid", "Fibroblasts", "Endothelial", "T cells", "B cells")
tmpfile <- paste0("plot/", dir.used, file.prefix, "PAS_VlnPlot_Imp.pdf")
tmp.text <- element_text(family="sans", size=16)
pdf(tmpfile, width = 12, height = 3.5)
path_plot_list <- lapply(dataset$DataSets[c(1,4,7:10)], function(tumor) {
    seu <- seu.list[[tumor]]
    DefaultAssay(seu) <- "Metab"
    seu <- seu[, seu$celltype %in% used_celltypes]
    features = sapply(used_path, grep, x = rownames(seu), fixed = T, value = T)
    p <- Seurat::VlnPlot(seu, features = features,
                          group.by = "celltype", cols = celltype_color, pt.size = 0) &
              geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
              plot_annotation(title = tumor) &
              theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, 
                    axis.line=element_line(size = .3, colour="black"))
    print(p)
    return(p)
})
dev.off()
# tmpfile <- metabte0("plot/", dir.used, file.prefix, "PAS_VlnPlot_Pooled_Imp.pdf")
# pdf(tmpfile, width = 12, height = 16)
# path_plot_list[[1]] <- 
# wrap_plots(tmp, nrow = 4)
# dev.off()


## 2. Metabolic genes --------
metab.gene <- split(gaudeMetabDf$GeneSymbol, gaudeMetabDf$Pathway)
used_gene <- c("GAPDH", metab.gene$`Glutamate metabolism`)
tmpfile <- paste0("plot/", dir.used, file.prefix, "Gene_VlnPlot.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 12, height = 20)
tmp <- lapply(dataset$DataSets[c(1,4,7:10)], function(tumor) {
    seu <- seu.list[[tumor]]
    DefaultAssay(seu) <- "RNA"
    print(Seurat::VlnPlot(seu, features = used_gene,
                          group.by = "celltype", cols = celltype_color, pt.size = 0) &
              geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
              plot_annotation(title = tumor) &
              theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text))
})
dev.off()


used_gene <- c("GLS", "GLUD1", "GLUL")
tmpfile <- paste0("plot/", dir.used, file.prefix, "Gene_VlnPlot_Imp.pdf")
tmp.text <- element_text(family="sans", size= 16)
pdf(tmpfile, width = 12, height = 3.5)
tmp <- lapply(dataset$DataSets[c(1,4,7:10)], function(tumor) {
    seu <- seu.list[[tumor]]
    seu <- seu[, seu$celltype %in% used_celltypes]
    DefaultAssay(seu) <- "RNA"
    p <- Seurat::VlnPlot(seu, features = used_gene,
                         group.by = "celltype", cols = celltype_color, pt.size = 0) &
        geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
        plot_annotation(title = tumor) &
        theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, 
              axis.line=element_line(size = .3, colour="black"))
    print(p)
    return(p)
})
dev.off()

# tmpfile <- metabte0("plot/", dir.used, file.prefix, "Gene_VlnPlot_Imp1.pdf")
# pdf(tmpfile, width = 12, height = 3.5)
# tmp[[1]] <- tmp[[1]] & xlab(NULL) & theme(axis.text.x = element_blank())
# tmp[[2]] <- tmp[[2]] & ggtitle(NULL) & xlab(NULL) & theme(axis.text.x = element_blank())
# tmp[[3]] <- tmp[[3]] & ggtitle(NULL) & xlab(NULL) & theme(axis.text.x = element_blank())
# tmp[[4]] <- tmp[[4]] & ggtitle(NULL) & xlab(NULL) 
# wrap_plots(tmp, ncol = 1)
# dev.off()
