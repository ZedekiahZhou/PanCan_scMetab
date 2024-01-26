# ==================================================================================================
# Author: Zhou Zhe
# Goals: Defferential analysis of 3 HDAC2 knockdown vs 3 WT H1299 cells lines
# Version: 1.0
# Date: Oct 28, 2022
# ==================================================================================================

rm(list=ls())
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggVolcano)
library(scMetab)
library(ggrepel)

# sample info
sam <- read.delim("src/05_epi/5.0_TF_RNAseq/config", header = FALSE)
colnames(sam) <- c("Sample", "Type", "Gene")
TERM2GENE <- data.frame(TERM = gaudeMetabDf$Pathway, GENE = gaudeMetabDf$GeneSymbol)
MP_list <- readRDS("res/02_NMF/MPs/MP_list.RDS")$Malignant
MP_df <- reshape2::melt(MP_list, measure.vars = colnames(MP_list))
colnames(MP_df) = c("TERM", "GENE")
TERM2GENE = MP_df


# # Prepare counts from featureCounts ----------------------------------------------------------------
# count1 <- read.delim("res/3_epi/RNAseq/counts/H1299-scramble-1_featurecounts.txt", 
#                      comment.char = "#", row.names = 1)
# count1 <- count1[, "gene_name", drop = FALSE]
# 
# # read count matrix
# for (i in 1:nrow(sam)) {
#     tmp <- read.delim(paste0("res/3_epi/RNAseq/counts/", sam$Sample[i], "_featurecounts.txt"),
#                       comment.char = "#", row.names = 1)
#     if (all(rownames(tmp) == rownames(count1))) {
#         count1[[sam$Sample[i]]] <- tmp[[7]]
#     }
# }
# write.table(count1, file = "res/3_epi/RNAseq/raw_counts_featureCounts.tsv", quote = F, sep = "\t")

# --------------------------------------------------------------------------------------------------
# DESeq2 ......
# --------------------------------------------------------------------------------------------------
# prepare for DESeq2
cts <- read.delim("res/05_epi/RNAseq/raw_counts_featureCounts.tsv", row.names = 1, check.names = FALSE)
id2symbol <- cts[, "gene_name", drop = FALSE]
cts <- cts[-1]

coldata <- data.frame(Type = factor(sam$Type, levels = c("WT", "KD")),
                      Gene = factor(sam$Gene, levels = c("HDAC2", "FOSL1")),
                      Sample = sam$Sample, row.names = sam$Sample)
hdac2 <- coldata$Sample[coldata$Gene == "HDAC2"]
fosl1 <- coldata$Sample[coldata$Gene == "FOSL1"]


## construct DESeq2 Obj for HDAC2 ------
DEfun <- function(samples, fout) {
    dds <- DESeqDataSetFromMatrix(countData = cts[, samples], colData = coldata[samples, ], design = ~ Type)
    dds <- dds[rowSums(counts(dds)) >= 10, ] # xxxx genes
    print(paste("There were", dim(dds)[1], "genes kept after genes with low expression were removed!"))
    
    dds <- DESeq(dds) # resultsNames(dds)
    res <- results(dds, name = "Type_KD_vs_WT", alpha = 0.05, independentFiltering = FALSE)
    print(summary(res))
    resOut <- data.frame(as.data.frame(res), symbol = id2symbol[rownames(res), ])
    
    rld <- assay(rlog(dds, blind = FALSE))
    resOut <- cbind(rld, resOut)
    write.csv(resOut, file = file.path("res/04_epi/RNAseq/", fout))
    return(resOut)
}


res_hdac2 <- DEfun(hdac2, "DEres_H1299_HDAC2_KO_v3.csv")
res_fosl1 <- DEfun(fosl1, "DEres_PANC1_FOSL1_KO_v3.csv")


GSEAfun <- function(res) {
    gl <- res$stat
    names(gl) <- res$symbol
    gl <- sort(gl, decreasing = TRUE)
    
    # go <- gseGO(gl, ont = "BP", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", eps = 0)
    metab  <- GSEA(gl, minGSSize = 5, TERM2GENE = TERM2GENE, pvalueCutoff = 1)
    return(list(metab = metab))
}


multiple_plot <- function(res, fout, gene, 
                          log2FC_cut = 0.4, padj_cut = 0.05) {
    
    tmp.text <- element_text(family="sans", size=6)
    
    ## sample dist heat map 
    sampleDists <- dist(t(res[1:6]))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- colnames(res)[1:6]
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pdf(file = paste0("plot/04_epi/RNAseq/", fout, "_dist_heatmap.pdf"), width = 6, height = 6)
    print(pheatmap(sampleDistMatrix, col=colors,
             clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists))
    dev.off()
    
    ## gsea plot
    gsea_res <- GSEAfun(res)
    metab_path <- data.frame(gsea_res$metab) %>% 
        filter(pvalue < 0.05) %>% arrange(NES) %>%
        mutate(ID = factor(.data$ID, levels = unique(.data$ID)), 
               Regulate = ifelse(NES > 0, "Up", "Down"))
    ggplot(metab_path, mapping = aes(x = ID, y = NES, fill = Regulate)) + 
        geom_bar(stat = "identity") + 
        scale_fill_brewer(palette = "Set1", direction = -1) + 
        theme_classic() + coord_flip() + 
        theme(text = tmp.text, axis.text.x = tmp.text, axis.text.y = tmp.text, legend.position = "top", 
              axis.line = element_line(linewidth = 0.5/.pt), axis.ticks = element_line(linewidth = 0.5/.pt), 
              legend.text = tmp.text, legend.key.size = unit(3, "mm"))
    ggsave(file = paste0("plot/04_epi/RNAseq/", fout, "_Hist.pdf"), width = 3, height = 1.5)
    
    ## Volcano plot
    Imp_genes <- unique(unlist(strsplit(metab_path$core_enrichment, "/", fixed = T)))
    Imp_genes <- data.frame(Gene = Imp_genes, pvalue = res$pvalue[match(Imp_genes, res$symbol)]) %>%
        arrange(pvalue)
    res <- res %>%
        # mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
        mutate(Regulate = ifelse(padj < padj_cut & log2FoldChange < -log2FC_cut, "Down", 
                                 ifelse(padj < padj_cut & log2FoldChange > log2FC_cut, "Up", "Normal")), 
               log10FDR = -log10(padj))
    Imp_genes <- Imp_genes %>% filter(Gene %in% res$symbol[res$Regulate != "Normal"])
    res <- res %>% mutate(label = ifelse(symbol %in% c(Imp_genes$Gene[1:10], gene), symbol, NA))
               
    mycol = brewer.pal(9, "Set1")[c(1, 2, 9)]
    names(mycol) = c("Up", "Down", "Normal")
    ggplot(res, mapping = aes(x = log2FoldChange, y = log10FDR, color = Regulate)) + 
        geom_point(size = 0.3) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.3) + 
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) + 
        geom_label_repel(aes(label = label), size = 6/.pt, max.overlaps = 1000, show.legend = F, seed = 200) + 
        scale_color_manual(values = mycol) + 
        xlab("log2(Fold Change)") + ylab("log10(FDR)") + 
        xlim(-5, 5) + theme_bw() + 
        theme(text = tmp.text, axis.text.x = tmp.text, axis.text.y = tmp.text, 
              legend.text = tmp.text)
    ggsave(file = paste0("plot/04_epi/RNAseq/", fout, "_Volcano.pdf"), width = 3.3, height = 2.5)
}

multiple_plot(res_fosl1, "Plot_PANC1_FOSL1_KO_v2", "FOSL1")
multiple_plot(res_hdac2, "Plot_H1299_HDAC2_KO_v2", "HDAC2")

