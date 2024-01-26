rm(list=ls())
library(data.table)
library(scMetab)
library(GSVA)
library(ggplot2)
library(ggrepel)

# I. Load data ==============
expr <- data.frame(fread("data/CCLE/CCLE_expression_full.csv"), row.names = 1, check.names = F)
tmp_gene <- colnames(expr)

#### update gene symbol
gene <- data.frame(tmpID = tmp_gene, symbol = sub(" \\(.+$", "", tmp_gene), ensembl = sub("\\)$", "", sub(".+ \\(", "", tmp_gene)))
gene$updated <- update_symbols(gene$symbol)$updated

#### subset genes
gene <- gene[!duplicated(gene$updated), ]
expr <- expr[, gene$tmpID]
colnames(expr) <- gene$updated
expr <- data.frame(t(expr), check.names = F)

#### read sample info
all_sample <- read.csv("data/CCLE/sample_info.csv")
brca_sample <- subset(all_sample, primary_disease == "Breast Cancer" & primary_or_metastasis == "Primary" &
                          lineage_sub_subtype %in% c("ERneg_HER2neg", "ERpos_HER2neg", "ERpos_HER2pos"))  # 22 cell lines
brca_sample$subtype = ifelse(brca_sample$lineage_sub_subtype == "ERneg_HER2neg", "TNBC", "ER_BC")
brca_sample <- brca_sample[, c("DepMap_ID", "CCLE_Name", "Subtype", "age", "lineage_sub_subtype", "lineage_molecular_subtype",
                               "patient_id", "Cellosaurus_NCIt_disease", "subtype")]

#### read metabolomics sample info
metab_sample <- readxl::read_excel("data/CCLE/SuppTables.xlsx", sheet = "1-cell line annotations")
all(brca_sample$CCLE_Name %in% metab_sample$`Name with tissue origins`)
brca_sample <- subset(brca_sample, CCLE_Name %in% metab_sample$`Name with tissue origins`)
#### remove 4 cell lines, finally 18 cell lines with both expression and metabolomics data

#### subset expr by sample
all(brca_sample$DepMap_ID %in% colnames(expr))
expr <- expr[, brca_sample$DepMap_ID]


#### gene level
res_gene <- plyr::ldply(data.frame(t(expr), check.names = FALSE), function(x) {
    tres = t.test(x[brca_sample$subtype == "ER_BC"], x[brca_sample$subtype == "TNBC"])
    data.frame(statistic = tres$statistic, pvalue = tres$p.value)
})
res_gene$padj <- p.adjust(res_gene$pvalue, method = "BH")



#### GSVA
metab.gene.list <- split(gaudeMetabDf$GeneSymbol, f = gaudeMetabDf$Pathway)
metab.gene.list <- metab.gene.list[sapply(metab.gene.list, length) >=5]

metab.gsva <- gsva(as.matrix(expr), metab.gene.list, method = "gsva", kcdf = "Gaussian", mx.diff = T, parallel.sz = 10)
metab.gsva <- data.frame(t(metab.gsva), check.names = F)
fastSave::saveRDS.pigz(metab.gsva, file = "data/CCLE/gene_expr_metab_gsva.rds")

brca_sample$DepMap_ID == rownames(metab.gsva)
res_gsva <- plyr::ldply(metab.gsva, function(x) {
    tres = t.test(x[brca_sample$subtype == "ER_BC"], x[brca_sample$subtype == "TNBC"])
    data.frame(statistic = tres$statistic, pvalue = tres$p.value)
})
res_gsva$padj <- p.adjust(res_gsva$pvalue, method = "BH")

# vertical
gsva_hist_data <- subset(res_gsva, padj < 0.05)
gsva_hist_data <- gsva_hist_data[order(gsva_hist_data$statistic), ]
gsva_hist_data$.id <- factor(gsva_hist_data$.id, levels = gsva_hist_data$.id)
gsva_hist_data$type = ifelse(gsva_hist_data$statistic > 0, "ERBC", "TNBC")

tumorType_color <- celltype_color[c(1, 5)]
names(tumorType_color) <- c("TNBC", "ERBC")

tmp.text <- element_text(family="sans", size=8)
p1 <- ggplot(gsva_hist_data, mapping = aes(x = .id, y = statistic, fill = type)) +
    geom_bar(stat = "identity", width = 0.6) +
    coord_flip() +
    # scale_fill_manual(values = brewer.pal(2, "RdYlBu")) +
    ggtitle("ER_BC vs TNBC") +
    xlab("Pathway") + ylab("t value") +
    scale_fill_manual(values = tumorType_color) +
    theme_classic() +
    theme(text = tmp.text, axis.text = tmp.text, legend.position = "none", plot.title = tmp.text,
          # axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line=element_line(size = .3, colour="black"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
ggsave(p1, filename = paste0("res/3_epi/CCLE/CCLE_BRCA_subtype_hist.pdf"), width = 5, height = 2)


#### read metabolomics data
metabolomics <- data.frame(readxl::read_excel("data/CCLE/SuppTables.xlsx", sheet = "1-clean data"), row.names = 1, check.names = F)
metabolomics <- data.frame(t(metabolomics), check.names = F)
metabolomics <- metabolomics[, brca_sample$CCLE_Name]


#### metabolomics
res_metab <- plyr::ldply(data.frame(t(metabolomics), check.names = FALSE), function(x) {
    tres = t.test(x[brca_sample$subtype == "ER_BC"], x[brca_sample$subtype == "TNBC"])
    d = mean(x[brca_sample$subtype == "ER_BC"]) - mean(x[brca_sample$subtype == "TNBC"])
    data.frame(statistic = tres$statistic, pvalue = tres$p.value, log2FC = log2(10^d))
})
res_metab$padj <- p.adjust(res_metab$pvalue, method = "BH")
res_metab$log10p <- -log10(res_metab$pvalue)
res_metab$sig <- ifelse(res_metab$pvalue >= 0.05, "Nonsignificant", ifelse(res_metab$statistic > 0, "ERBC", "TNBC"))
res_metab$ID <- ifelse(res_metab$pvalue < 0.05, res_metab$.id, "")
write.csv(res_metab, file = "res/3_epi/CCLE/ER_v_TNBC_metabolites.csv", row.names = F, quote = F)

res_metab$log10p <- ifelse(res_metab$log10p > 4, 4, res_metab$log10p)
res_metab$log2FC <- ifelse(res_metab$log2FC < -4, -4, res_metab$log2FC)
tmp.text <- element_text(family="sans", size=8)
tumorType_color <- c(celltype_color[c(1, 5)], "grey")
names(tumorType_color) <- c("TNBC", "ERBC", "Nonsignificant")


p2 = ggplot(res_metab, aes(x = log2FC, y = log10p, color = sig, label = ID)) +
    geom_point() +
    geom_text_repel(size = 2.88, max.overlaps = 100, show.legend = F) +
    xlab("log2 FC") + ylab("-log10(P value)") +
    xlim(c(-4, 4)) + ylim(c(0, 4)) +
    scale_color_manual(values = tumorType_color) +
    # geom_vline(xintercept=-log10(0.05)) + geom_hline(yintercept=-log10(0.05)) +
    theme_classic() +
    theme(axis.text=tmp.text, text = tmp.text, legend.position = "bottom")
ggsave(p2, filename = paste0("res/3_epi/CCLE/BRCA_metabolites_volcano_changed_color.pdf"), width = 4, height = 4)
