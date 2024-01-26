# ====================================================================
# Author: Zhou Zhe
# Function: Dot plot for metabolic pathway changes of
#               specific celltypes (Epithelial, myeloidblast, Myeloid) between tumor and normal
#               include cancer: TNBC, LUAD, COAD, PAAD
# Version: 2.0
# Date: Jul 18, 2022
# ====================================================================

rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(patchwork)
library(tidyverse)
library(scMetab)
library(ggsci)

pseudo <- 1e-50

dir.used <- "2_TvN/DE/"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
gaude.list <- split(gaudeMetabDf$GeneSymbol, f = gaudeMetabDf$Pathway)

# All Cell Types =====================================================
## 1. Gene markers ----
celltypes <- c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")
gene_markers.list <- lapply(celltypes, function(used_celltype) {
    used_celltype = sub(" ", "_", used_celltype)
    print(paste(Sys.time(), "--", used_celltype))
    file.prefix <- paste0("Celltype_TvN_", used_celltype, "_")
    
    gene_markers <- lapply(c(1, 3:10), function(i) {
        tumor = dataset$DataSets[i]
        tmp.data <- read.delim(paste0("res/2_TvN/AUCell/", tumor, "/", tumor, "_", used_celltype, "_Gene_Markers.tsv"))
        tmp.data$tumor <- tumor
        tmp.data
    })
    
    gene_markers <- do.call(rbind, gene_markers)
    gene_markers$cluster = used_celltype
    gene_markers <- subset(gene_markers, (p_val_adj < 0.01) & (abs(avg_log2FC) > 0.1) & ((pct.1 > 0.1) | (pct.2 > 0.1)))
    gene_markers$log10padj <- (-1) * sign(gene_markers$avg_log2FC) * log10(gene_markers$p_val_adj)
    gene_markers$log10padj_pseudo <- (-1) * sign(gene_markers$avg_log2FC) * log10(gene_markers$p_val_adj + pseudo)
    
    gene_markers <- gene_markers[order(gene_markers$avg_log2FC, decreasing = T), ]
    gene_markers$gene <- factor(gene_markers$gene, levels = unique(gene_markers$gene))
    
    top10 <- gene_markers %>%
        group_by(tumor) %>%
        top_n(n = -10, wt = p_val_adj)
    top10 <- top10[order(top10$tumor, top10$p_val_adj), ]
    top10$gene <- factor(top10$gene, levels = unique(top10$gene))
    
    tmpfile <- paste0("plot/", dir.used, file.prefix, "Gene_dotplot.pdf")
    tmp.text <- element_text(family="sans", size=8)
    mycolor <- brewer.pal(11, "RdBu")[c(10, 6, 2)]
    pdf(tmpfile, width = 3.5, height = ceiling(length(unique(top10$gene))/9))
    print(ggplot(top10, aes(x = tumor, y = gene, color = log10padj_pseudo)) + 
              geom_point(aes(size = abs(avg_log2FC))) + 
              scale_size(limits = c(0, max(abs(top10$avg_log2FC))), range = c(0.5, 3.5)) + 
              scale_color_gradient2(low = mycolor[1], high = mycolor[3], mid = mycolor[2]) + 
              theme(axis.text=tmp.text, text = tmp.text, 
                    axis.text.x = element_text(angle = 45, hjust = 1), 
                    axis.line=element_line(size = .3, colour="black")))
    dev.off()
    return(gene_markers)
})
names(gene_markers.list) <- sub(" ", "_", celltypes)
TvN_gene_markers <- do.call(rbind, gene_markers.list)
write.table(TvN_gene_markers, file = paste0("res/", dir.used, "TvN_Gene_Markers_All_Datasets.tsv"), 
            quote = F, sep = "\t", row.names = F)
TvN_metab_gene_markers <- subset(TvN_gene_markers, gene %in% gaudeMetabDf$GeneSymbol)
write.table(TvN_metab_gene_markers, file = paste0("res/", dir.used, "TvN_Metab_Gene_Markers_All_Datasets.tsv"), 
            quote = F, sep = "\t", row.names = F)

### Number of Up-Down datasets ------
TvN_MetabGene_markers <- subset(TvN_gene_markers, gene %in% gaudeMetabDf$GeneSymbol)
gene_marker_list <- split(TvN_MetabGene_markers, f = TvN_MetabGene_markers$cluster)
updown_list <- lapply(gene_marker_list, function(x) {
    x <- x %>%
        group_by(gene) %>%
        summarise(up = sum(avg_log2FC > 0), down = sum(avg_log2FC < 0), 
                  up_m_down = sum(avg_log2FC > 0) - sum(avg_log2FC < 0), 
                  up_p_down = n(), label = "")
})
updown_list$Malignant$label = ifelse(updown_list$Malignant$up_p_down > 5, 
                                     updown_list$Malignant$gene, "")
updown_list$Myeloid$label = ifelse(updown_list$Myeloid$up_p_down > 5, 
                                   updown_list$Myeloid$gene, "")
updown_list$T_cells$label = ifelse(updown_list$T_cells$up_p_down >= 4, 
                                   updown_list$T_cells$gene, "")
updown_list$B_cells$label = ifelse(updown_list$B_cells$up_p_down >= 4, 
                                   updown_list$B_cells$gene, "")
updown_list$Endothelial$label = ifelse(updown_list$Endothelial$up_p_down >= 7, 
                                   updown_list$Endothelial$gene, "")
updown_list$Fibroblasts$label = ifelse(updown_list$Fibroblasts$up_p_down >= 7, 
                                   updown_list$Fibroblasts$gene, "")

tmp.text <- element_text(family="sans", size=8)
pdf(paste0("plot/", dir.used, "TvN_Gene_UpDown_counts_with_label.pdf"), width = 20, height = 8)
for (celltype in names(updown_list)) {
    plot_data <- updown_list[[celltype]]
    p = ggplot(plot_data, aes(up_p_down, up_m_down, color = up_m_down)) + 
        geom_count() + ylab("N(up) - N(down)") + xlab("N(up) + N(down)") + 
        ggrepel::geom_text_repel(aes(up_p_down, up_m_down, label = label), 
                                 size = 8*5/14, max.overlaps = 10) + 
        scale_color_gradient2(high = brewer.pal(11, "RdYlBu")[2], low = brewer.pal(11, "RdYlBu")[10], 
                              mid = "grey", midpoint = 0) +
        theme_classic() + ggtitle(celltype) + 
        theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, 
              plot.title = tmp.text,
              axis.line=element_line(size = .3, colour="black")) 
    print(p)
}
dev.off()


## 2. Metabolic markers ----
celltypes <- c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")
metab_markers.list <- lapply(celltypes, function(used_celltype) {
    used_celltype = sub(" ", "_", used_celltype)
    print(paste(Sys.time(), "--", used_celltype))
    file.prefix <- paste0("Celltype_TvN_", used_celltype, "_")
    
    metab_markers <- lapply(c(1, 3:10), function(i) {
        tumor = dataset$DataSets[i]
        tmp.data <- read.delim(paste0("res/2_TvN/AUCell/", tumor, "/", tumor, "_", used_celltype, "_Metab_Markers.tsv"))
        tmp.data$tumor <- tumor
        tmp.data
    })
    metab_markers <- do.call(rbind, metab_markers)
    metab_markers$pathway <- gaude_trans$updated[match(metab_markers$pathway, 
                                                       gaude_trans$original)]
    metab_markers$cluster <- used_celltype
    metab_markers <- subset(metab_markers, (p_val_adj < 0.01) & (abs(avg_log2FC) > 0.01) & ((pct.1 > 0.1) | (pct.2 > 0.1)))
    
    metab_markers$log10padj <- sign(metab_markers$avg_log2FC) * (-1) * log10(metab_markers$p_val_adj)
    metab_markers$log10padj_pseudo <- sign(metab_markers$avg_log2FC) * (-1) * log10(metab_markers$p_val_adj + pseudo)
    metab_markers <- metab_markers[order(metab_markers$avg_log2FC, decreasing = T), ]
    metab_markers$pathway <- factor(metab_markers$pathway, levels = unique(metab_markers$pathway))
    ### core enrichment
    gene_markers <- gene_markers.list[[used_celltype]]
    metab_markers$core <- apply(metab_markers[, c("tumor", "pathway")], 1, function(x) {
        genes = intersect(gaude.list[[x[2]]], gene_markers$gene[gene_markers$tumor == x[1]])
        return(paste0(genes, collapse = "/"))
    })
    
    tmpfile <- paste0("plot/", dir.used, file.prefix, "Metab_dotplot.pdf")
    tmp.text <- element_text(family="sans", size=8)
    pdf(tmpfile, width = 5.5, height = ceiling(length(unique(metab_markers$pathway))/9))
    print(ggplot(metab_markers, aes(x = tumor, y = pathway, color = log10padj_pseudo)) + 
        geom_point(aes(size = abs(avg_log2FC))) + 
        scale_color_distiller(palette = "RdBu") + 
        scale_size(limits = c(0, 0.5), range = c(0.5, 3.5)) + 
        theme_classic() + 
        theme(axis.text=tmp.text, text = tmp.text, 
              axis.text.x = element_text(angle = 45, hjust = 1), 
              axis.line=element_line(size = .3, colour="black")))
    dev.off()
    return(metab_markers)
})
names(metab_markers.list) <- sub(" ", "_", celltypes)
TvN_metab_markers <- do.call(rbind, metab_markers.list)
TvN_metab_markers$pathway <- gaude_trans$updated[match(TvN_metab_markers$pathway, 
                                                  gaude_trans$original)]
write.table(TvN_metab_markers, file = paste0("res/", dir.used, "TvN_Metab_Markers_All_Datasets.tsv"), 
            quote = F, sep = "\t", row.names = F)


## 3. TF markers ----
celltypes <- c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")
TF_markers.list <- lapply(celltypes, function(used_celltype) {
    used_celltype = sub(" ", "_", used_celltype)
    print(paste(Sys.time(), "--", used_celltype))
    file.prefix <- paste0("Celltype_TvN_", used_celltype, "_")
    
    TF_markers <- lapply(c(1, 3:10), function(i) {
        tumor = dataset$DataSets[i]
        tmp.data <- read.delim(paste0("res/2_TvN/AUCell/", tumor, "/", tumor, "_", used_celltype, "_TF_Markers.tsv"))
        tmp.data$tumor <- tumor
        tmp.data
    })
    TF_markers <- do.call(rbind, TF_markers)
    TF_markers$cluster <- used_celltype
    TF_markers <- subset(TF_markers, (p_val_adj < 0.01) & (abs(avg_log2FC) > 0.01) & ((pct.1 > 0.1) | (pct.2 > 0.1)))
    
    TF_markers$log10padj <- sign(TF_markers$avg_log2FC) * (-1) * log10(TF_markers$p_val_adj)
    TF_markers$log10padj_pseudo <- sign(TF_markers$avg_log2FC) * (-1) * log10(TF_markers$p_val_adj + pseudo)
    TF_markers <- TF_markers[order(TF_markers$avg_log2FC, decreasing = T), ]
    TF_markers$pathway <- factor(TF_markers$pathway, levels = unique(TF_markers$pathway))
    ### core enrichment
    gene_markers <- gene_markers.list[[used_celltype]]
    
    regulons <- lapply(c(1, 3:10), function(i) {
        tumor = dataset$DataSets[i]
        geneIds(getGmt(paste0("res/5_SCENIC/regulon/", tumor, "_reg_metab.gmt")))
    })
    names(regulons) = dataset$DataSets[c(1, 3:10)]
    
    TF_markers$core <- apply(TF_markers[, c("tumor", "pathway")], 1, function(x) {
        genes = intersect(regulons[[x[1]]][[x[2]]], gene_markers$gene[gene_markers$tumor == x[1]])
        return(paste0(genes, collapse = "/"))
    })
    
    tmpfile <- paste0("plot/", dir.used, file.prefix, "TF_dotplot.pdf")
    tmp.text <- element_text(family="sans", size=8)
    pdf(tmpfile, width = 3.5, height = ceiling(length(unique(TF_markers$pathway))/9))
    print(ggplot(TF_markers, aes(x = tumor, y = pathway, color = log10padj_pseudo)) + 
              geom_point(aes(size = abs(avg_log2FC))) + 
              scale_color_distiller(palette = "RdBu") + 
              scale_size(limits = c(0, 0.5), range = c(0.5, 3.5)) + 
              theme_classic() + 
              theme(axis.text=tmp.text, text = tmp.text, 
                    axis.text.x = element_text(angle = 45, hjust = 1), 
                    axis.line=element_line(size = .3, colour="black")))
    dev.off()
    return(TF_markers)
})
names(TF_markers.list) <- celltypes
TvN_TF_markers <- do.call(rbind, TF_markers.list)
write.table(TvN_TF_markers, file = paste0("res/", dir.used, "TvN_TF_Markers_All_Datasets.tsv"), 
            quote = F, sep = "\t", row.names = F)




# All combined ====================================================================================
TvN_metab_markers <- read.delim("res/2_TvN/DE/TvN_Metab_Markers_All_Datasets.tsv")
TvN_metab_markers <- TvN_metab_markers[order(TvN_metab_markers$avg_log2FC, decreasing = T), ]
TvN_metab_markers$pathway <- factor(TvN_metab_markers$pathway, levels = unique(TvN_metab_markers$pathway))
tmp.text <- element_text(family="sans", size=8)

tmpfile <- paste0("plot/", dir.used, "ALL_pooled_dotplot.pdf")
pdf(tmpfile, width = 12, height = length(levels(TvN_metab_markers$pathway))/9)
print(ggplot(TvN_metab_markers, aes(x = tumor, y = pathway, color = log10padj_pseudo)) + 
          geom_point(aes(size = abs(avg_log2FC))) + 
          scale_color_distiller(palette = "RdBu") + 
          scale_size(limits = c(0, 0.5), range = c(0.5, 3.5)) + 
          theme_classic() + 
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black")) + 
          facet_wrap(~ cluster, nrow = 1))
dev.off()

# Immune cells ==============
immune <- subset(TvN_metab_markers, cluster %in% c("Myeloid", "T_cells", "B_cells"))
immune$pathway <- as.character(immune$pathway)
immune <- immune[order(immune$avg_log2FC, decreasing = T), ]
immune$pathway <- factor(immune$pathway, levels = unique(immune$pathway))

tmpfile <- paste0("plot/", dir.used, "Immune_pooled_dotplot.pdf")
pdf(tmpfile, width = 8, height = length(levels(immune$pathway))/9)
print(ggplot(immune, aes(x = tumor, y = pathway, color = log10padj_pseudo)) + 
          geom_point(aes(size = abs(avg_log2FC))) + 
          scale_color_distiller(palette = "RdBu") + 
          scale_size(limits = c(0, 0.5), range = c(0.5, 3.5)) + 
          theme_classic() + 
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black")) + 
          facet_wrap(~ cluster, nrow = 1))
dev.off()


# Strommal cells ==============
strommal <- subset(TvN_metab_markers, cluster %in% c("Endothelial", "Fibroblasts"))
strommal$pathway <- as.character(strommal$pathway)
strommal <- strommal[order(strommal$avg_log2FC, decreasing = T), ]
strommal$pathway <- factor(strommal$pathway, levels = unique(strommal$pathway))

tmpfile <- paste0("plot/", dir.used, "Strommal_pooled_dotplot.pdf")
pdf(tmpfile, width = 7, height = length(levels(strommal$pathway))/9)
print(ggplot(strommal, aes(x = tumor, y = pathway, color = log10padj_pseudo)) + 
          geom_point(aes(size = abs(avg_log2FC))) + 
          scale_color_distiller(palette = "RdBu") + 
          scale_size(limits = c(0, 0.5), range = c(0.5, 3.5)) + 
          theme_classic() + 
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black")) + 
          facet_wrap(~ cluster, nrow = 1))
dev.off()

tmpfile <- paste0("plot/", dir.used, "Strommal_dotplot.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 6, height = length(levels(strommal$pathway))/9)
print(ggplot(strommal, aes(x = celltype_tumor, y = pathway, color = log10padj)) + 
          geom_point(aes(size = abs(avg_diff))) + 
          scale_color_distiller(palette = "RdBu") + 
          scale_size(limits = c(0, 0.5), range = c(0.5, 3.5)) + 
          theme_classic() + 
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black")))
print(ggplot(strommal, aes(x = celltype_tumor, y = pathway, color = log10padj)) + 
          geom_point(aes(size = abs(avg_log2FC))) + 
          scale_color_distiller(palette = "RdBu") + 
          scale_size(limits = c(0, 0.5), range = c(0.5, 3.5)) + 
          theme_classic() + 
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black")))
dev.off()