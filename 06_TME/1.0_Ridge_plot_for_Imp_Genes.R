rm(list=ls())
library(dplyr)
library(Seurat)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggridges)

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets


# CD4 T cells
MP_list = readRDS("res/02_NMF/MPs/MP_list.RDS")$CD4_T
Imp_path <- c("CD4.c20(TNFRSF9+ Treg)")
Imp_gene <- c(MP_list$MMP3, "COX5A", "GAPDH")

sig_score_per_sample <- readRDS(file = paste0("res/02_NMF/MPs/CD4_T_sig_score_per_sample.RDS"))
sig_score_per_sample <- lapply(sig_score_per_sample, scale)
sig_score = do.call(rbind, sig_score_per_sample)

expr_per_sample = readRDS("data/rds/NMF/log_transformed/CD4_T_log2_all_samples.RDS")
expr = do.call(rbind, lapply(expr_per_sample, function(x) t(as.matrix(x[Imp_gene, ]))))

df = data.frame(cbind(expr, sig_score), check.names = F)

res.list <- lapply(colnames(expr), function(x) {
    df$group = ifelse(df[[x]] > 0, paste0(x, "+"), paste0(x, "-"))
    ggplot(df, aes(y = group, x = .data[[Imp_path]], fill = group)) + 
        geom_density_ridges() + 
        scale_fill_manual(values = brewer.pal(9, "Set1")[2:1]) + 
        theme_bw() + 
        theme(panel.grid = element_blank())
})
wrap_plots(res.list, ncol = 4) + plot_annotation("CD4_T")
ggsave("plot/06_T/CD4T_cor/CD4_Ridge_Imp_gene.pdf", width = 16, height = 12)


# CD8 T cells
MP_list = readRDS("res/02_NMF/MPs/MP_list.RDS")$CD8_T
Imp_path <- c("CD8.c12(terminal Tex)")
Imp_gene <- c(MP_list$MMP4, "COX5A", "GAPDH")

sig_score_per_sample <- readRDS(file = paste0("res/02_NMF/MPs/CD8_T_sig_score_per_sample.RDS"))
sig_score_per_sample <- lapply(sig_score_per_sample, scale)
sig_score = do.call(rbind, sig_score_per_sample)

expr_per_sample = readRDS("data/rds/NMF/log_transformed/CD8_T_log2_all_samples.RDS")
expr = do.call(rbind, lapply(expr_per_sample, function(x) t(as.matrix(x[Imp_gene, ]))))

df = data.frame(cbind(expr, sig_score), check.names = F)

res.list <- lapply(colnames(expr), function(x) {
    df$group = ifelse(df[[x]] > 0, paste0(x, "+"), paste0(x, "-"))
    ggplot(df, aes(y = group, x = .data[[Imp_path]], fill = group)) + 
        geom_density_ridges() + 
        scale_fill_manual(values = brewer.pal(9, "Set1")[2:1]) + 
        theme_bw() + 
        theme(panel.grid = element_blank())
})
wrap_plots(res.list, ncol = 4) + plot_annotation("CD8_T")
ggsave("plot/06_T/CD8T_cor/CD8_Ridge_Imp_gene.pdf", width = 16, height = 12)

# Myeloid cells
rm(list=ls())
MP_list = readRDS("res/02_NMF/MPs/MP_list.RDS")$Myeloid
Imp_path <- c("TIM_M2", "TIM_Phagocytosis")
Imp_gene <- c(MP_list$MMP1, "COX5A", "GAPDH", "GLUL")

sig_score_per_sample <- readRDS(file = paste0("res/02_NMF/MPs/Myeloid_sig_score_per_sample.RDS"))
sig_score_per_sample <- lapply(sig_score_per_sample, scale)
sig_score = do.call(rbind, sig_score_per_sample)

expr_per_sample = readRDS("data/rds/NMF/log_transformed/Myeloid_log2_all_samples.RDS")
expr = do.call(rbind, lapply(expr_per_sample, function(x) t(as.matrix(x[Imp_gene, ]))))

df = data.frame(cbind(expr, sig_score), check.names = F)

pdf("plot/06_T/Mye/Mye_Ridge_Imp_gene.pdf", width = 16, height = 13.5)
for (i in 1:2) {
    res.list <- lapply(colnames(expr), function(x) {
        df$group = ifelse(df[[x]] > 0, paste0(x, "+"), paste0(x, "-"))
        ggplot(df, aes(y = group, x = .data[[Imp_path[i]]], fill = group)) + 
            geom_density_ridges() + 
            scale_fill_manual(values = brewer.pal(9, "Set1")[2:1]) + 
            theme_bw() + 
            theme(panel.grid = element_blank())
    })
    print(wrap_plots(res.list, ncol = 4) + plot_annotation("Myeloid"))
}
dev.off()
