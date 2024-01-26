library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Expression vs Variability ===========================

## Expression
Malignant_cpm_all_samples <- readRDS("data/rds/NMF/cpm/Malignant_cpm_all_samples.RDS")
MP_list <- readRDS("res/02_NMF/MPs/MP_list.RDS")$Malignant
sam_info_list <- readRDS("data/rds/NMF/cells_info_per_sam.RDS")$Malignant
attach(sam_info_list)

MP_expr <- sapply(Malignant_cpm_all_samples, function(x) {
    tmp <- sapply(MP_list, function(MP) {
        MP <- MP[MP %in% rownames(x)]
        colMeans(as.matrix(x[MP, ]))
    })
    tmp <- log2(tmp/10 + 1)
    tmp <- colMeans(tmp)
})
MP_expr <- data.frame(t(MP_expr), check.names = F)
all(rownames(MP_expr) == sam_info$samID)
MP_expr <- do.call(rbind, tapply(MP_expr, sam_info$Dataset, colMeans))


## Variability
load("res/02_NMF/MPs/Malignant_MPs_p30_.7_.2_final.RData")
Cluster_sams <- lapply(Cluster_list, function(x) unique(sub("_rank4_9.+", "", x)))
sam_list <- split(sam_info$samID, sam_info$Dataset)
MP_var <- sapply(Cluster_sams, function(x) {
    sapply(sam_list, function(sams) {
        sum(x %in% sams)/length(sams)
    })
})
colnames(MP_var) <- colnames(MP_expr)


## plot
plot(MP_expr[, 1], MP_var[, 1])

tmp_text = element_text(family = "sans", size = 6)
plist = lapply(1:15, function(i) {
    df = data.frame(expr = MP_expr[, i], varibility = MP_var[, i], tumor = rownames(MP_expr))
    p = ggplot(df, aes(x = expr, y = varibility)) + 
        geom_point(size = 0.3) + 
        ggpubr::stat_cor(method = "pearson", size = 6/.pt) +
        geom_smooth(method = lm, se = F, linetype = "dashed", color = "red4", linewidth = 0.3) +
        theme_classic() + 
        geom_text_repel(aes(label = tumor), size = 6/.pt) + 
        ggtitle(label = paste0("MMP", i)) + 
        theme(plot.title = element_text(family = "sans", size = 6, face = "bold"),
            axis.title = element_blank(), axis.text = tmp_text, axis.line = element_line(linewidth = 0.3))
})
p = wrap_plots(plist, ncol = 4)
ggsave("plot/02_NMF/Malignant/Malignant_MP_expr_vs_varibility_Comb.pdf", p, width = 6, height = 6)

pdf("plot/02_NMF/Malignant_MP_expr_vs_varibility.pdf")

dev.off()



# 画MP9 expression heatmap
# 看MMP9的样本间差异