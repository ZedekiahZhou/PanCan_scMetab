# ====================================================================
# Author: Zhou Zhe
# Function: Inter-sample correlation with clinical info
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

dir.used <- "04_epi/InterSample/"
dir.create(file.path("res", dir.used), recursive = T)
dir.create(file.path("plot", dir.used), recursive = T)
# file.prefix <- "patient_celltype"

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
seed.use = 2021


# I. calculate CNV score and ITH ==================================================================
samCNV = list()
for (i in 1:10) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    seu_obj <- readRDS(paste0("data/rds/MetabGene/", tumor, ".TpN.unintegrated.rds"))
    if (!("samID" %in% colnames(seu_obj[[]]))) seu_obj$samID = seu_obj$patientID
    if (tumor == "PRAD") seu_obj$samID <- seu_obj$samID2
     c 
    if (tumor == "BRCA_valid2") {
        infercnv <- readRDS("../other2/EMBO_2020_Wu_TNBC/res/0_reproduce/infercnv_server/run.final.infercnv_obj")
        cnv_score1 <- colSums((infercnv@expr.data - 1)^2)
        infercnv <- readRDS("../other2/WangShu_BRCA//res/0_reproduce/infercnv_server/run.final.infercnv_obj")
        cnv_score2 <- colSums((infercnv@expr.data - 1)^2)
        cnv_score <- c(cnv_score1, cnv_score2)
    } else {
        infercnv <- readRDS(paste0("../", dataset$Directory[i],"/res/0_reproduce/infercnv_server/run.final.infercnv_obj"))
        cnv_score <- colSums((infercnv@expr.data - 1)^2)
    }
    
    seu_obj$cnv_score <- cnv_score
    
    epi <- subset(seu_obj, celltype == "Malignant" & TorN == "T")
    
    pca.data <- data.frame(Embeddings(epi, reduction = "pca")[, 1:30])
    all(rownames(pca.data) == colnames(epi))
    pca.list <- split(pca.data, epi$samID)
    ITH <- lapply(pca.list, function(x) {
        m = colMeans(x)
        delta = apply(x, 2, sd)
        y = sweep(x, 2, m)
        
        outlier = sweep(abs(y), 2, 3*delta) > 0
        idx <- rowSums(outlier[, 1:3]) == 3
        
        res = rowSums(y^2)
        res[idx] = NA
        return(res)
    })
    names(ITH) <- NULL
    ITH <- unlist(ITH)
    all(colnames(epi) == names(ITH))
    epi$ITH <- ITH[match(colnames(epi), names(ITH))]
    
    saminfo <- epi[[c("samID", "cnv_score", "ITH")]]
    saminfo <- saminfo %>%
        group_by(samID) %>%
        summarise(cnv_score = mean(cnv_score, na.rm = T), ITH = mean(ITH, na.rm = T))
    samCNV[[tumor]] <- saminfo
}
samCNV <- do.call(rbind, samCNV)
write.table(samCNV, file = "res/3_epi/InterSample/sam_CNV_ITH.tsv", sep = "\t", row.names = F)


# II. Pathway correaltion with clinical info ======================================================
pooled.merged <- readRDS(file = "data/rds/pseudo_sam_celltype.rds")
pooled.merged <- subset(pooled.merged, celltype == "Malignant" & n >= 20)
metab.gsva <- readRDS("data/rds/pseudo_sam_celltype_metab_gsva.rds")
metab.gsva <- metab.gsva[, colnames(pooled.merged)]
pooled.merged[["GSVA"]] <- CreateAssayObject(counts = metab.gsva)

sam_clin <- read.delim("data/Pooled_clinical_info.tsv", na.strings = c("NA", ""))
sam_clin <- merge(sam_clin, samCNV, by = "samID", all = T)
sam_clin <- sam_clin[match(pooled.merged$samID, sam_clin$samID), ]
rownames(sam_clin) <- paste(sam_clin$dataset, sam_clin$samID, "Malignant", sep = "_")

pathway_cor <- lapply(1:10, function(i) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    seu_obj <-pooled.merged[, pooled.merged$tumor == tumor]
    gsva <- as.data.frame(t(as.matrix(GetAssayData(seu_obj, assay = "GSVA", slot = "counts"))))
    saminfo <- sam_clin[rownames(gsva), ]
    
    clin_cor <- lapply(c("n", "subtype", "age", "gender", "stage", "cnv_score", "ITH"), function(var) {
        if (length(unique(saminfo[[var]])) > 1) {
            res <- plyr::ldply(gsva, function(path) {
                res.lm <- lm(path~saminfo[[var]])
                res.lm.s <- summary(res.lm)
                data.frame(dataset = tumor,
                           cate=var, 
                           r.squared = res.lm.s$r.squared, 
                           adj.r.squared = res.lm.s$adj.r.squared, 
                           p.value = anova(res.lm)$'Pr(>F)'[1])
            }) 
            res$p.adj <- p.adjust(res$p.value, method = "BH")
            return(res)
        } else {
            return(data.frame(.id = colnames(gsva), 
                              dataset = tumor,
                              cate = var, 
                              r.squared = NA, 
                              adj.r.squared = NA,
                              p.value = NA, 
                              p.adj = NA))
        }
        
    })
    clin_cor = do.call(rbind, clin_cor)
    # names(clin_cor) <- c("n", "subtype", "age", "gender", "stage")
    return(clin_cor)
})
# names(pathway_cor) <- dataset$DataSets
pathway_cor = do.call(rbind, pathway_cor)
pathway_cor$cate <- ifelse(pathway_cor$cate == "cnv_score", "cnvScore", pathway_cor$cate)
write.table(pathway_cor, file = "res/3_epi/InterSample/sample_clinical_cor.tsv", quote = F, sep = "\t", row.names = F)
# saveRDS(pathway_cor, file = "res/3_epi/InterSample/sample_clinical_cor.rds")


# III. Plot =======================================================================================
pathway_cor <- read.delim("res/04_epi/InterSample/sample_clinical_cor.tsv")
pathway_cor$.id <- gaude_trans[pathway_cor$.id, ]$updated
# remove LUAD_valid1 and variable n
# pathway_cor <- subset(pathway_cor, cate != "n" & !(dataset %in% c("LUAD_valid1", "PRAD")))
pathway_cor <- subset(pathway_cor, cate != "n" & !(dataset %in% c("LUAD_valid1")))
pathway_cor <- subset(pathway_cor, cate != "subtype")

r2adj <- dcast(pathway_cor, .id~dataset+cate, value.var = "adj.r.squared") %>% column_to_rownames(".id")
pvalue <- dcast(pathway_cor, .id~dataset+cate, value.var = "p.value") %>% column_to_rownames(".id")
padj <- dcast(pathway_cor, .id~dataset+cate, value.var = "p.adj") %>% column_to_rownames(".id")


## remove ITH columns
r2adj <- r2adj[, grep("ITH", colnames(r2adj), value = T, invert = T)]
pvalue <- pvalue[, grep("ITH", colnames(pvalue), value = T, invert = T)]

# keep correlation > 0.1 and adjust p < 0.05
plot_data1 = plot_data2 = r2adj
idx1 <- (is.na(pvalue)) | (pvalue > 0.05) | (is.na(r2adj)) | (abs(r2adj) < 0.1)
plot_data1[idx1] <- 0
# plot_data1 <- plot_data1[rowSums(plot_data1 != 0) > 4, ]

nSigs <- sapply(split.default(plot_data1, rep(1:4, 9)), function(x) rowSums(x != 0))
idx <- rowSums(nSigs > 3) > 0
plot_data1 = plot_data1[idx, ]


column_annot <- sub("_[^_]+$", "", colnames(r2adj))
myCols <- rev(brewer.pal(n=5, name="RdBu"))[c(3, 5)]
tmp.text <- gpar(fontfamily="sans", fontsize = 6)
pdf(paste0("plot/", dir.used, "PanCan_intersample_clinical_heatmap_change_color.pdf"), width = 8, height = 8)
Heatmap(plot_data1, col = circlize::colorRamp2(c(0, 1), myCols), 
        rect_gp = gpar(col = "white"), 
        cluster_columns = FALSE, cluster_rows = TRUE, 
        row_names_side = "left", show_row_dend = FALSE, 
        row_names_max_width = unit(10, "cm"),
        row_names_gp = tmp.text, column_names_gp = tmp.text, 
        row_title_gp = tmp.text, column_title_gp = tmp.text, 
        column_title_rot = 45, height = unit(nrow(plot_data1)/9, "inches"), 
        width = unit(1.5*ncol(plot_data1)/10, "cm"), 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        top_annotation = HeatmapAnnotation(Tumor = column_annot, col = list(Tumor = dataset_color), 
                                           simple_anno_size = unit(2, "mm"), 
                                           annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
                                           which = "column"), 
        column_split = sub("^.+_", "", colnames(r2adj)), gap = unit(1, "mm"), 
        show_column_names = F, column_names_side = "top", 
        column_names_max_height = unit(6, "cm"),
        border = TRUE)
dev.off()




## Pathway with tumor type (Only for test)
pathway_cor <- read.delim("res/3_epi/InterSample/sample_clinical_cor.tsv")
# remove LUAD_valid1, PRAD and variable n
pathway_cor <- subset(pathway_cor, cate == "subtype" & !(dataset %in% c("LUAD_valid1", "PRAD")))

r2adj <- dcast(pathway_cor, .id~dataset+cate, value.var = "adj.r.squared") %>% column_to_rownames(".id")
pvalue <- dcast(pathway_cor, .id~dataset+cate, value.var = "p.value") %>% column_to_rownames(".id")
padj <- dcast(pathway_cor, .id~dataset+cate, value.var = "p.adj") %>% column_to_rownames(".id")

# keep correlation > 0.1 and adjust p < 0.05
plot_data1 = plot_data2 = r2adj
idx1 <- (is.na(pvalue)) | (pvalue > 0.05) | (is.na(r2adj)) | (abs(r2adj) < 0.1)
plot_data1[idx1] <- 0
# plot_data1 <- plot_data1[rowSums(plot_data1 != 0) > 0, ]

nSigs <- sapply(split.default(plot_data1, sub("_.+", "", colnames(cor_df_comnbined))), 
                function(x) rowSums(x != 0))
idx <- rowSums(nSigs > 3) > 0
plot_data = cor_df_comnbined[idx, ]

column_annot <- sub("_[^_]+$", "", colnames(r2adj))
myCols <- rev(brewer.pal(n=5, name="RdBu"))[c(3, 5)]
tmp.text <- gpar(fontfamily="sans", fontsize = 8)
pdf(paste0("plot/", dir.used, "PanCan_intersample_cancerSubType_heatmap.pdf"), width = 6, height = 6)
Heatmap(plot_data1, col = circlize::colorRamp2(c(0, 1), myCols), 
        rect_gp = gpar(col = "white"), 
        cluster_columns = FALSE, cluster_rows = TRUE, 
        row_names_side = "left", show_row_dend = FALSE, 
        row_names_max_width = unit(10, "cm"),
        row_names_gp = tmp.text, column_names_gp = tmp.text, 
        row_title_gp = tmp.text, column_title_gp = tmp.text, 
        heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
        top_annotation = HeatmapAnnotation(Tumor = column_annot, col = list(Tumor = dataset_color), 
                                           simple_anno_size = unit(2, "mm"), 
                                           annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
                                           which = "column"), 
        column_split = sub("^.+_", "", colnames(r2adj)), gap = unit(1, "mm"), 
        show_column_names = F, column_names_side = "top", 
        column_names_max_height = unit(6, "cm"),
        border = TRUE)
dev.off()
