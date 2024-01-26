# prostate cancer have only 9 cell lines (<10), so removed. 
rm(list=ls())

library(scMetab)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(scatterpie)
library(GSVA)
library(RColorBrewer)
library(dplyr)

# I. Prepare data ========
## read drug information
drug <- readxl::read_excel("data/CCLE/Cell_Iorio_2016_TableS1_Drug_Info.xlsx", 
                           sheet = "TableS1F_ScreenedCompounds", skip = 2)
drug$drug_ID = paste0("ID_", drug$Identifier)
drug = drug[drug$`Clinical Stage` != "experimental", ] ## only keep drugs "clinically approved" and "in clinical development"

## read IC50 data
Drug2ID = data.frame(t(readxl::read_excel("data/CCLE/Cell_Iorio_2016_TableS4_CCLE_drug_IC50.xlsx", 
                                             sheet = 1, skip = 4, na = "NA", n_max = 1)[, c(-1, -2)]))
colnames(Drug2ID) = "Name"
Drug2ID$drug_ID = paste0("ID_", rownames(Drug2ID))
IC50 <- data.frame(readxl::read_excel("data/CCLE/Cell_Iorio_2016_TableS4_CCLE_drug_IC50.xlsx", 
                           sheet = 1, skip = 5, na = "NA"), check.names = F)
colnames(IC50) = c("COSMICID", "cl_name", paste0("ID_", rownames(Drug2ID)))

Drug2ID = Drug2ID[match(drug$drug_ID, Drug2ID$drug_ID), ]
IC50 = IC50[, c("COSMICID", "cl_name", drug$drug_ID)]

## read gene expression
expr <- data.frame(data.table::fread("data/CCLE/CCLE_expression_full.csv"), row.names = 1, check.names = F)
all_sample <- read.csv("data/CCLE/sample_info.csv") # 1840 cell lines
kept_cl = all_sample[(all_sample$COSMICID %in% IC50$COSMICID) & (all_sample$DepMap_ID %in% rownames(expr)), ] # 695 cell lines left

## update gene symbol
gene <- data.frame(tmpID = colnames(expr), symbol = sub(" \\(.+$", "", colnames(expr)), 
                   ensembl = sub("\\)$", "", sub(".+ \\(", "", colnames(expr))))
gene$updated <- update_symbols(gene$symbol)$updated

## remove genes with duplicated names
gene <- gene[!duplicated(gene$updated), ]
expr <- expr[, gene$tmpID]
colnames(expr) <- gene$updated

## subset expr matrix for used cell lines
sub_expr <- expr[kept_cl$DepMap_ID, ]
rownames(sub_expr) <- kept_cl$stripped_cell_line_name
sub_expr <- data.frame(t(sub_expr), check.names = F)

## subset IC50 data 
IC50 <- IC50[match(kept_cl$COSMICID, IC50$COSMICID), c(-1,-2)]
rownames(IC50) <- kept_cl$stripped_cell_line_name

# II. GSVA ========
MP_list = as.list(readRDS("res/02_NMF/MPs/MP_list.RDS")$Malignant)
# MP_gsva <- gsva(as.matrix(sub_expr), MP_list, method = "gsva", kcdf = "Gaussian", mx.diff = T, parallel.sz = 10)
# MP_gsva <- data.frame(t(MP_gsva), check.names = F)
# save(MP_gsva, kept_cl, file = "res/02_NMF/MP_gsva.rda")
load("res/02_NMF/MP_gsva.rda")

## correlation
cl_list = split(kept_cl, kept_cl$primary_disease) 
cl_list <- cl_list[sapply(cl_list, nrow)>=10]
cor_pearson_list = lapply(cl_list, function(x) {
    Hmisc::rcorr(as.matrix(MP_gsva[x$stripped_cell_line_name, ]), 
                 as.matrix(IC50[x$stripped_cell_line_name, ]), type = "pearson")
})
cor_spearman_list = lapply(cl_list, function(x) {
    Hmisc::rcorr(as.matrix(MP_gsva[x$stripped_cell_line_name, ]), 
                 as.matrix(IC50[x$stripped_cell_line_name, ]), type = "spearman")
})

## prepare data for plot
melt_fun <- function(x, slot, value.name) {
    tmp_df = data.frame(x[[slot]][grep("MMP", rownames(x[[slot]]), invert = T), 
                                  grep("MMP", colnames(x[[slot]]))], check.names = F)
    tmp_df$Drugs = rownames(tmp_df)
    res1 = reshape2::melt(tmp_df, measure.vars = setdiff(colnames(tmp_df), "Drugs"), id.vars = "Drugs",
                          variable.name = "MMPs", value.name = value.name)
}
r_pearson = lapply(cor_pearson_list, melt_fun, slot = "r", value.name = "R_pearson")
p_pearson = lapply(cor_pearson_list, melt_fun, slot = "P", value.name = "P_pearson")
r_spearman = lapply(cor_spearman_list, melt_fun, slot = "r", value.name = "R_spearman")
p_spearman = lapply(cor_spearman_list, melt_fun, slot = "P", value.name = "P_spearman")
cor_df = lapply(1:length(r_pearson), function(i) {
    df = Reduce(function(x, y) dplyr::left_join(x, y, by = c("Drugs", "MMPs")), 
                list(r_pearson[[i]], p_pearson[[i]], r_spearman[[i]], p_spearman[[i]]))
    df$Padj_pearson = p.adjust(df$P_pearson, method = "BH")
    df$Padj_spearman = p.adjust(df$P_spearman, method = "BH")
    return(df)
})
names(cor_df) = names(r_pearson)

# III. Plot ==================
### scatter pie plot
imp_df = cor_df[c("Breast Cancer", "Colon/Colorectal Cancer", "Gastric Cancer", "Lung Cancer", "Pancreatic Cancer")]
comb_data = plyr::ldply(names(imp_df), function(tumor) {
    x = imp_df[[tumor]]
    x$Cancer = tumor
    x$drug_name = Drug2ID$Name[match(x$Drugs, Drug2ID$drug_ID)]
    return(x)
})
comb_data$value = 1
sig_comb_data <- filter(comb_data, !is.na(R_pearson) & Padj_pearson < 0.05 & Padj_spearman < 0.05 & abs(R_pearson) > 0.3)

used_drugs = unique(sig_comb_data$Drugs)
used_MMPs = unique(as.character(sig_comb_data$MMPs))
sig_comb_data$Drugs = match(sig_comb_data$Drugs, used_drugs)
sig_comb_data$MMPs = match(sig_comb_data$MMPs, used_MMPs)

tmp.text <- element_text(family="sans", size=6)
ggplot() + 
    geom_scatterpie(aes(x=Drugs, y=MMPs), data=sig_comb_data, cols="Cancer", long_format=TRUE, size = 0.3, pie_scale = 0.8) + 
    theme_bw() + 
    scale_x_continuous(breaks = seq(1, length(used_drugs), 1), expand = c(0.01, 0.01), labels = Drug2ID$Name[match(used_drugs, Drug2ID$drug_ID)]) + 
    scale_y_continuous(breaks = seq(1, length(used_MMPs), 1), expand = c(0.01, 0.01), labels = used_MMPs) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, axis.title = tmp.text) + coord_fixed()
ggplot() + 
    geom_scatterpie(aes(x=Drugs, y=MMPs, r = -log(P_pearson)/40), data=sig_comb_data, cols="Cancer", long_format=TRUE, 
                    size = 0.3, pie_scale = 0.8) + 
    theme_bw() + 
    scale_x_continuous(breaks = seq(1, length(used_drugs), 1), expand = c(0.01, 0.01), labels = Drug2ID$Name[match(used_drugs, Drug2ID$drug_ID)]) + 
    scale_y_continuous(breaks = seq(1, length(used_MMPs), 1), expand = c(0.01, 0.01), labels = used_MMPs) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, axis.title = tmp.text) + coord_fixed()

ggsave(filename = "plot/02_NMF/Drugs/Drugs_scatterpie.pdf", width = 6, height = 4)


# IV. Elesclomol (ID_1031) volvano ======
comb_data = dplyr::mutate(comb_data, 
                          Significant = ifelse(!is.na(R_pearson) & Padj_pearson < 0.05 & Padj_spearman < 0.05 & abs(R_pearson) > 0.3, 
                                               "Significant", "Non-Significant"))
comb_data = dplyr::mutate(comb_data, Group = ifelse(Significant == "Significant" & drug_name == "Elesclomol", 
                                                    ifelse(R_pearson > 0, "positive", "negative"), "None"))
comb_data$Label = ifelse(comb_data$Group != "None", paste0(comb_data$MMPs, " - ", comb_data$Cancer), "")
tmp_color = c(celltype_color[1:2], "black")
names(tmp_color) = c("positive", "negative", "None")
ggplot(comb_data, aes(x = R_pearson, y = -log10(Padj_pearson), color = Group, label = Label)) + 
    geom_point() + 
    scale_color_manual(values = tmp_color) + 
    geom_hline(yintercept = -log10(0.05), color = "red4", linetype = "dashed") + 
    xlab("R") + ylab("-log10(FDR)") + 
    geom_text_repel(max.overlaps = 100) + 
    theme_classic()
ggsave("plot/02_NMF/Drugs/All_drug_valcono.pdf", width = 4, height = 4)

tmp <- lapply(cor_df, function(x) x[x$Drugs == "ID_1031", ])
a = reshape2::melt(tmp, id.vars = colnames(tmp[[1]]))
a = dplyr::mutate(a, 
                Significant = ifelse(!is.na(R_pearson) & Padj_pearson < 0.05 & Padj_spearman < 0.05 & abs(R_pearson) > 0.3, 
                                     "Significant", "Non-Significant"))
a$MMPs = sub(" .+", "", a$MMPs)
a$tissue = sub(".+/", "", sub(" Cancer", "", a$L1))
a = dplyr::mutate(a, Group = ifelse(Significant == "Significant" , ifelse(R_pearson > 0, "positive", "negative"), "None"))
a$Label = ifelse(a$Group != "None", paste0(sub(" Cancer", "", a$tissue), " - ", a$MMPs), "")
# a = comb_data[comb_data$drug_name == "Elesclomol", ]

tmp_color = c(celltype_color[1:2], "black")
names(tmp_color) = c("positive", "negative", "None")
ggplot(a, aes(x = R_pearson, y = -log10(Padj_pearson), color = Group, label = Label)) + 
    geom_hline(yintercept = -log10(0.05), color = "red4", linetype = "dashed", linewidth = 0.3) + 
    geom_vline(xintercept = -0.3, color = "red4", linetype = "dashed", linewidth = 0.3) + 
    geom_vline(xintercept = 0.3, color = "red4", linetype = "dashed", linewidth = 0.3) + 
    geom_point(size = 0.3) + 
    scale_color_manual(values = tmp_color) + 
    scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) + 
    xlab("R") + ylab("-log10(FDR)") + 
    geom_text_repel(size = 6/.pt) + 
    theme_classic() + 
    theme(text=tmp.text, legend.position = "None", 
          axis.line = element_line(linewidth = 0.3), axis.text = tmp.text, axis.title = tmp.text)
ggsave("plot/02_NMF/Drugs/Elescolomol_only_Valcono.pdf", width = 2, height = 2.5)


# V. Elesclomol IC50 vs MMP7 score ======
cl_coad = cl_list$`Colon/Colorectal Cancer`$stripped_cell_line_name
plot(MP_gsva[cl_coad, ]$`MMP7 Glycosylation1`, IC50[cl_coad, ]$ID_1031)
cor.test(MP_gsva[cl_coad, ]$`MMP7 Glycosylation1`, IC50[cl_coad, ]$ID_1031)

tmp1 <- data.frame(MMP7_score = MP_gsva[cl_coad, ]$`MMP7 Glycosylation1`, 
                   IC50_Elesclomol = IC50[cl_coad, ]$ID_1031)
ggplot(tmp1, aes(x = MMP7_score, y = IC50_Elesclomol)) + 
    geom_point(size = 0.3) + 
    geom_smooth(method = lm, se = F, linetype = "dashed", color = "red4", linewidth = 0.3) +
    xlab("MMP7 Score") + ylab("IC50 (Elesclomol)") + 
    stat_cor(method = "pearson", size = 6/.pt) + 
    scale_x_continuous(breaks = c(-0.5, 0, 0.5)) + 
    theme_classic() + 
    theme(text=tmp.text, legend.position = "None", 
          axis.line = element_line(linewidth = 0.3), axis.text = tmp.text, axis.title = tmp.text)
ggsave("plot/02_NMF/Drugs/CRC_ELesclomol_IC50_vs_MMP7.pdf", width = 1.5, height = 1.25)


# V. Axitinib IC50 vs MMP4 score ======
cl_luad = cl_list$`Lung Cancer`$stripped_cell_line_name
plot(MP_gsva[cl_luad, ]$`MMP4 Glycolysis`, IC50[cl_luad, ]$ID_1021)
cor.test(MP_gsva[cl_luad, ]$`MMP4 Glycolysis`, IC50[cl_luad, ]$ID_1021)

tmp2 <- data.frame(MMP4_score = MP_gsva[cl_luad, ]$`MMP4 Glycolysis`, 
                   IC50_Axitinib = IC50[cl_luad, ]$ID_1021)
ggplot(tmp2[tmp2$IC50_Axitinib > 0, ], aes(x = MMP4_score, y = IC50_Axitinib)) + 
    geom_point(size = 0.3) + 
    geom_smooth(method = lm, se = F, linetype = "dashed", color = "red4", linewidth = 0.3) +
    xlab("MMP4 Score") + ylab("IC50 (Axitinib)") + 
    stat_cor(method = "pearson", size = 6/.pt) + 
    scale_x_continuous(breaks = c(-0.5, 0, 0.5)) + 
    theme_classic() + 
    theme(text=tmp.text, legend.position = "None", 
          axis.line = element_line(linewidth = 0.3), axis.text = tmp.text, axis.title = tmp.text)
ggsave("plot/02_NMF/Drugs/LUAD_Axitinib_IC50_vs_MMP4.pdf", width = 1.5, height = 1.25)





# VI. Experiment validation ========
read_dat <- function(Gene) {
    dat <- read.delim(paste0("res/02_NMF/Elescolomol/", Gene, ".txt"), skip = 3, header = F, nrows = 1)[, c(-1, -2, -99)]
    dat <- matrix(unlist(dat), nrow = 8, byrow = T)[2:7, 2:8]
    dat_norm <- data.frame(dat/mean(dat[, 1]))
    colnames(dat_norm) <- c(0, as.character(10^c(-3:2)))
    dat_df = reshape2::melt(dat_norm)
    dat_df$Group = Gene
    #dat_df$variable = factor(dat_df$variable, levels = c("Inf", as.character(-3:2)))
    return(dat_df)
}

used_gene = "FUT3+AOC1+TST"
treat = read_dat(used_gene)
empty = read_dat("Ctrl")
df = rbind(treat, empty)
df = df[df$variable %in% c("0", "0.01", "1", "100"), ]
# a$variable = as.numeric(as.character(a$variable))

tmp_color = c(celltype_color[1:2])
names(tmp_color) = c(used_gene, "Ctrl")
tmp.text <- element_text(family="sans", size=6)
p = ggplot(df, aes(x = variable, y = value, group = Group, color = Group)) + 
    stat_summary(fun = "mean", geom = "point", size = 0.1) +
    stat_summary(fun = "mean", geom = "line", linewidth = 0.3) + 
    scale_color_manual(values = tmp_color) + 
    # scale_x_discrete(breaks = c("Inf", "-2", "0", "2"), labels = c("0", "0.01", "0", "100")) + 
    # geom_smooth(method = NULL, se = F) +  
    ylab("Fold Viability") + xlab("Elesclomol (uM)") + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, linewidth = 0.3) + 
    stat_compare_means(method = "t.test", show.legend = F, label = "p.signif", size = 6/.pt) + 
    theme_classic() + 
    theme(text=tmp.text, legend.position = "top", legend.margin = margin(0.5,0.5,0.5,0.5), 
          legend.justification = "top",
          legend.box.spacing = unit(1, "pt"), legend.text = tmp.text, plot.margin = margin(0.5,0.5,0.5,0.5),
          axis.line = element_line(linewidth = 0.3), axis.text = tmp.text, axis.title = tmp.text)
ggsave(paste0("plot/02_NMF/Drugs/Elescolomol_HCT116_", used_gene, ".pdf"), p, width = 1.5, height = 1.25)


# 添加星号
# 确认用t还是wilcox














### heatmap plot for test
my_cols = brewer.pal(11, "Spectral")[c(10, 6, 2)]
tmp.text <- element_text(family="sans", size=8)
plot_fun <- function(M) {
    M = M[!is.na(M$R_pearson) & M$Padj_pearson < 0.05 & M$Padj_spearman < 0.05 & abs(M$R_pearson) > 0.5, ]
    M = M[order(M$P_pearson), ]
    M$Drugs = factor(M$Drugs, levels = unique(M$Drugs))
    ggplot(M, aes(x = Drugs, y = MMPs, fill = R_pearson)) + 
        geom_tile() + 
        scale_color_gradient2(low = my_cols[1], mid = my_cols[2], high = my_cols[3], na.value = my_cols[2]) + 
        scale_fill_gradient2(low = my_cols[1], mid = my_cols[2], high = my_cols[3], na.value = my_cols[2]) + 
        theme_bw() +
        scale_x_discrete(expand = c(0,0)) + 
        scale_y_discrete(expand = c(0,0)) + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), 
              axis.ticks = element_blank(), text = tmp.text, 
              legend.key = element_rect(colour = "black", linewidth = 0.5), 
              legend.key.size = unit(4, "mm"))
    
}

plot_fun(M = cor_df$`Pancreatic Cancer`)
plot_fun(M = cor_df$`Colon/Colorectal Cancer`)
plot_fun(M = cor_df$`Gastric Cancer`)




### important genes for MMP7 drug sensitivity
df_coad <- cbind(MP_gsva[cl_coad, , drop = F], 
                 IC50[cl_coad, "Elesclomol", drop = F], 
                 t(sub_expr)[cl_coad, MP_list$`MMP7 Glycosylation1`])            
a = Hmisc::rcorr(as.matrix(df_coad))

plot(df_coad$`MMP7 Glycosylation1`, df_coad$CKB)
plot(df_coad$`MMP7 Glycosylation1`, df_coad$TST)



