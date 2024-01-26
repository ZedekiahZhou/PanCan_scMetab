# ==================================================================================================
# Author: Zhou Zhe
# Program: TCGA-GTEx data TvN using wilcox test
# Version: 1.0
# Date: Feb 14, 2022
# ==================================================================================================

###### ====== Prepare TCGA data ===== ######
rm(list=ls())
library(clusterProfiler)
library(Seurat)
library(scMetab)
library(ggsci)
library(ggplot2)
library(parallel)
cl = makeCluster(10)

dir.used <- "2_TvN/TCGA/"
dir.create(file.path("plot/", dir.used), recursive = T)
dir.create(file.path("res/", dir.used), recursive = T)
mycolor <- RColorBrewer::brewer.pal(11, "Spectral")[c(2,10)]

# I. Prepare data ========
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
dataset <- unique(dataset[, c("Tissue", "tumorType")])

# TCGA sample phenotype
pheno <- read.delim(gzfile(paste0("data/xena_tpm/TcgaTargetGTEX_phenotype.txt.gz")))
pheno$patient <- sub("-[0-9]+$", "", pheno$sample)
tissues = c("Breast", "Colon", "Rectum", "Lung", "Pancreas", "Prostate", "Stomach")
pheno <- subset(pheno, X_primary_site %in% tissues)
pheno$tissue <- ifelse(pheno$X_primary_site %in% c("Colon", "Rectum"), "Colorectum", pheno$X_primary_site)
pheno$TorN <- ifelse(pheno$X_sample_type == "Primary Tumor", "T", "N")

table(pheno$X_sample_type, pheno$X_study)
pheno <- pheno[pheno$X_sample_type %in% c("Normal Tissue", "Primary Tumor", "Solid Tissue Normal"), ]

# ID to Symbol
id2symbol <- read.delim("data/xena_tpm/probeMap_gencode.v23.annotation.gene.probemap")
id2symbol$updated <- update_symbols(id2symbol$gene)$updated
id2symbol <- id2symbol[!duplicated(id2symbol$updated), ]

# TCGA TvN
all_gsea <- mclapply(1:6, function(i) {
    tumor = dataset$tumorType[i]
    tpm_data <- data.frame(fread(paste0("data/xena_tpm/", dataset$Tissue[i], "_matrix.tsv.gz")), 
                           row.names = 1, check.names = F)
    
    tumor_pheno <- subset(pheno, tissue == dataset$Tissue[i])
    if (tumor == "LUAD") {
        tumor_pheno <- subset(pheno, detailed_category != "Lung Squamous Cell Carcinoma")
    }
    
    # only keep samples with both expression data and clinical data ---
    tumor_sam <- intersect(colnames(tpm_data), tumor_pheno$sample)
    print(paste("There are", ncol(tpm_data), "samples have expression data and", nrow(tumor_pheno), "samples have clinical info!"))
    print(paste("Among them there are", length(tumor_sam), "samples in common!"))
    tpm_data <- tpm_data[, tumor_sam]
    tumor_pheno <- tumor_pheno[match(tumor_sam, tumor_pheno$sample), ]
    samplesT <- tumor_pheno$sample[tumor_pheno$TorN == "T"]
    samplesN <- tumor_pheno$sample[tumor_pheno$TorN == "N"]
    
    # filter gene --
    tpm_data <- tpm_data[rownames(tpm_data) %in% id2symbol$id, ]
    rownames(tpm_data) <- id2symbol$updated[match(rownames(tpm_data), id2symbol$id)]
    
    # save file for TvN analysis
    fastSave::save.pigz(tpm_data, tumor_pheno, file = paste0("res/2_TvN/TCGA/", tumor, "_data_for_TvN.rda"), n.cores = 10)
    
    ###### ====== DE ====== ######
    tpm_data <- 2^tpm_data - 0.001
    res <- data.frame(Tumor_mean = rowMeans(tpm_data[, samplesT]),
                      Tumor_median = apply(tpm_data[, samplesT], 1, median), 
                      Normal_mean = rowMeans(tpm_data[, samplesN]), 
                      Normal_median = apply(tpm_data[, samplesN], 1, median))
    res$log2FC <- log2(res$Tumor_mean/res$Normal_mean)
    # clusterExport(cl, varlist = list("samplesT", "samplesN"))
    # system.time(res$p.value <- parApply(cl, tpm_data, 1, function(x) {wilcox.test(x[samplesT], x[samplesN])$p.value}))
    system.time(res$p.value <- apply(tpm_data, 1, function(x) {wilcox.test(x[samplesT], x[samplesN])$p.value}))
    res$stats <- log2(res$p.value) * (-1) * sign(res$log2FC)
    res$symbol <- rownames(res)
    
    ###### ====== GSEA ====== ######
    # Metabolic gene based GSEA
    TERM2GENE <- data.frame(TERM = gaudeMetabDf$Pathway, GENE = gaudeMetabDf$GeneSymbol)
    gl <- res$stats
    names(gl) <- rownames(res)
    gl <- gl[!is.na(gl)]
    
    ## Correction for promiscuity?? (not applied)
    # tmp <- as.data.frame(table(gaudeMetabDf$GeneSymbol))
    # promiscuity <- ifelse(names(gl) %in% tmp$Var1, yes = tmp$Freq, no = 1)
    # gl <- gl / promiscuity
    
    gl <- sort(gl, decreasing = TRUE)
    # gl[abs(gl) == Inf] <- sign(gl[abs(gl) == Inf]) * (max(abs(gl[abs(gl) != Inf])) + 1)
    summary(gl)
    gsea <- GSEA(gl, minGSSize = 5, TERM2GENE = TERM2GENE, pvalueCutoff = 1)
    gsea.df <- data.frame(gsea)
    
    
    ###### ====== Plot ====== ######
    tmpfile <- paste0("plot/", dir.used, tumor, "_GSEA_plot.pdf")
    tmp.text <- element_text(family="sans", size=17)
    pdf(tmpfile, width = 6, height = 5)
    for (pathway in gsea.df$ID) {
        print(enrichplot::gseaplot2(gsea, geneSetID = pathway, title = pathway) +
                  theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text) +
                  scale_fill_manual(values = mycolor[2:1]))
    }
    dev.off()

    write.table(res, file = paste0("res/", dir.used, tumor, "_Gene_markers.tsv"), quote = F, row.names = F, sep = "\t")
    write.table(gsea.df, file = paste0("res/", dir.used, tumor, "_PAS_GSEA.tsv"), quote = F, row.names = F, sep = "\t")
    
    return(gsea.df)
}, mc.cores = 6)


####### ================ TCGA TvN ==================== ###################
pas_markers <- lapply(1:6, function(i) {
    tumor <- dataset$tumorType[i]
    tmp.data <- read.delim(paste0("res/2_TvN/TCGA/", tumor, "_PAS_GSEA.tsv"))
    tmp.data$cancerType <- tumor
    tmp.data
})
pas_markers <- do.call(rbind, pas_markers)

pas_markers <- pas_markers[pas_markers$p.adjust < 0.05, ]

pas_markers$log10padj <- sign(pas_markers$NES) * (-1) * log10(pas_markers$p.adjust + 1e-3)
levels <- unique(pas_markers$ID[order(pas_markers$NES, decreasing = T)])
pas_markers$ID <- factor(pas_markers$ID, levels = levels)

tmpfile <- paste0("plot/", dir.used, "PanCan_dotplot.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(tmpfile, width = 5, height = 2.5)
print(ggplot(pas_markers, aes(x = cancerType, y = ID, color = log10padj)) + 
          geom_point(aes(size = abs(NES))) + 
          scale_color_distiller(palette = "RdBu") + 
          scale_size(limits = c(1, 3), range = c(1, 3.5)) + 
          theme_classic() +
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black")))
dev.off()