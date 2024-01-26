# ==================================================================================================
# Author: Zhou Zhe
# Program: Celltype TvN
# Version: 1.0
# Date: Jul 13, 2022
# ==================================================================================================

# init -------
rm(list=ls())

library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(ggplot2)
library(ggsci)
library(scMetab)
library(RColorBrewer)
# detach("package:scMetab", unload = TRUE)

n_cores <- 8
seed.use = 2021
TvN_color <- brewer.pal(11, "Spectral")[c(2,10)]

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
dataset <- dataset[-2, ]

## construct compare data.frame
celltypes <- c("B cells", "Endothelial", "Fibroblasts", "Myeloid", "T cells", "Mast cells")
compare.df <- data.frame(Tumor = paste0(celltypes, "_T"), Normal = paste0(celltypes, "_N"))
compare.df <- rbind(compare.df, data.frame(Tumor = "Malignant_T", Normal = "Epithelial_N"))
rownames(compare.df) <- c(celltypes, "Malignant")

for (i in 7:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
 
    # I. Load Data ========
    seu_obj <- readRDS(file = paste0("data/rds/PAS/", tumor, ".TpN.PAS.unintegrated.rds"))
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj <- NormalizeData(seu_obj)
    metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu_obj@assays$RNA))
    uniq_celltype = unique(seu_obj$celltype)
    
    # specify celltype used
    dir.used <- paste0("2_TvN/AUCell/", tumor, "/")
    dir.create(paste0("plot/", dir.used))
    dir.create(paste0("res/", dir.used), recursive = T)
    
    for (used_celltype in rownames(compare.df)) {
        if (!(used_celltype %in% uniq_celltype)) next
        
        print(paste(Sys.time(), "--", used_celltype))
        file.prefix <- paste0(tumor, "_", sub(" ", "_", used_celltype), "_")
        
        seu_celltype <- subset(seu_obj, celltype_x_TorN %in% unlist(compare.df[used_celltype, ]))
        table(seu_celltype$celltype_x_TorN)
        
       
        # II. DE ========
        DefaultAssay(seu_celltype)
        Idents(seu_celltype) <- "celltype_x_TorN"
        
        # metab markers
        metab_markers <- FindMarkers(seu_celltype, slot = "data", assay = "Metab",
                                   ident.1 = compare.df[used_celltype, "Tumor"],
                                   ident.2 = compare.df[used_celltype, "Normal"],
                                   logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                                   max.cells.per.ident = 2000, random.seed = seed.use)
        metab_markers$pathway <- rownames(metab_markers)
        
        # TF markers
        TF_markers <- FindMarkers(seu_celltype, slot = "data", assay = "TF",
                                     ident.1 = compare.df[used_celltype, "Tumor"],
                                     ident.2 = compare.df[used_celltype, "Normal"],
                                     logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                                     max.cells.per.ident = 2000, random.seed = seed.use)
        TF_markers$pathway <- rownames(TF_markers)
        
        
        # Metabolic gene markers
        gene_markers <- FindMarkers(seu_celltype, slot = "data", assay = "RNA", features = rownames(seu_celltype[["RNA"]]),
                                    ident.1 = compare.df[used_celltype, "Tumor"],
                                    ident.2 = compare.df[used_celltype, "Normal"],
                                    logfc.threshold = 0, min.pct = 0, only.pos = FALSE,
                                    max.cells.per.ident = 1000, random.seed = seed.use)
        gene_markers$gene <- rownames(gene_markers)
        
        
        # Metabolic gene based GSEA
        TERM2GENE <- data.frame(TERM = gaudeMetabDf$Pathway, GENE = gaudeMetabDf$GeneSymbol)
        gl <- log2(gene_markers$p_val) * (-1) * sign(gene_markers$avg_log2FC)
        names(gl) <- rownames(gene_markers)
        gl <- sort(gl, decreasing = TRUE)
        gl[abs(gl) == Inf] <- sign(gl[abs(gl) == Inf]) * (max(abs(gl[abs(gl) != Inf])) + 1)
        gsea <- GSEA(gl, minGSSize = 5, TERM2GENE = TERM2GENE, pvalueCutoff = 1, nPermSimple = 10000)
        gsea.df <- data.frame(gsea)
        
        tmpfile <- paste0("plot/", dir.used, file.prefix, "GSEA_plot.pdf")
        tmp.text <- element_text(family="sans", size=17)
        pdf(tmpfile, width = 6, height = 5)
        for (pathway in gsea.df$ID) {
            print(enrichplot::gseaplot2(gsea, geneSetID = pathway, title = pathway) +
                      theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text) +
                      scale_fill_manual(values = TvN_color[2:1]))
        }
        dev.off()
        
        
        tmpfile <- paste0("res/", dir.used, file.prefix)
        write.table(metab_markers, file = paste0(tmpfile, "Metab_Markers.tsv"), quote = F, row.names = F, sep = "\t")
        write.table(TF_markers, file = paste0(tmpfile, "TF_Markers.tsv"), quote = F, row.names = F, sep = "\t")
        write.table(gsea.df, file = paste0(tmpfile, "GSEA.tsv"), quote = F, row.names = F, sep = "\t")
        write.table(gene_markers, file = paste0(tmpfile, "Gene_Markers.tsv"), quote = F, row.names = F, sep = "\t")
    }
}
