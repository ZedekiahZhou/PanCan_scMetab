# ==================================================================================================
# Author: Zhou Zhe
# Program: Calculate and plot MP correlation with Signature path of different cell types
# Version: 1.0
# Date: Aug 13, 2023
# ==================================================================================================

rm(list=ls())
library(ComplexHeatmap)
library(scMetab)

# Prepare scores ===================================================================================
sig_path_list  <- list(Malignant = c("EMT", "Inflammatory response", "Angiogenesis", "Apoptosis", 
                                     "mTORC1 signaling", "DNA repair", "G2M checkpoint"), 
                       Myeloid = c("TIM_M1", "TIM_M2", "TIM_Angiogenesis", "TIM_Phagocytosis"), 
                       CD4_T = grep("CD4", unique(tcellSig$Pathway), value = T), 
                       CD8_T = grep("CD8", unique(tcellSig$Pathway), value = T))
celltype = "Malignant"
sig_path = sig_path_list[[celltype]]

#### read MP scores
MP_scores_per_sample <- readRDS(paste0("res/02_NMF/MPs/", celltype, "_MP_scores_p30_.7_.2_final_MinScore_0.8.RDS"))
sam_info_list <- readRDS("data/rds/NMF/cells_info_per_sam.RDS")[[celltype]]
attach(sam_info_list)

Gavish_scores_per_sample <- readRDS("res/02_NMF/MPs/Malignant_Gavish_MP_scores.RDS")
sig_score_per_sample <- readRDS(file = paste0("res/02_NMF/MPs/", celltype, "_sig_score_per_sample.RDS"))
TFreg_score_per_sample <- readRDS(file = paste0("res/02_NMF/MPs/", celltype, "_TFreg_score_per_sample.RDS"))
path_score_per_sample <- readRDS(file = paste0("res/02_NMF/MPs/", celltype, "_path_score_per_sample.RDS"))

# #### arrange sig scores
# sig_score <- lapply(unique(sam_info$Dataset), function(tumor) {
#     sig_pas <- data.frame(data.table::fread(paste0("res/02_SCENIC/aucell/", tumor, "_celltypeSig_aucell.csv")),
#                           row.names = 1, check.names = F)
#     sig_pas <- sig_pas[, sig_path]
#     
#     used_sams = sam_info$samID[sam_info$Dataset == tumor]
#     tmp_res <- lapply(used_sams, function(sam) {
#         sig_pas[sam_cells[[sam]], ]
#     })
#     names(tmp_res) <- used_sams
#     return(tmp_res)
# })
# sig_score_per_sample = do.call(c, sig_score)
# saveRDS(sig_score_per_sample, file = paste0("res/02_NMF/MPs/", celltype, "_sig_score_per_sample.RDS"))
# 
# #### arrange TF_reg scores
# TFreg_score <- lapply(unique(sam_info$Dataset), function(tumor) {
#     TFreg_pas <- data.frame(data.table::fread(paste0("res/02_SCENIC/aucell/", tumor, "_TF_aucell.csv")), 
#                             row.names = 1, check.names = F)
#     
#     used_sams = sam_info$samID[sam_info$Dataset == tumor]
#     tmp_res <- lapply(used_sams, function(sam) {
#         TFreg_pas[sam_cells[[sam]], ]
#     })
#     names(tmp_res) <- used_sams
#     return(tmp_res)
# })
# TFreg_score_per_sample <- do.call(c, TFreg_score)
# saveRDS(TFreg_score_per_sample, file = paste0("res/02_NMF/MPs/", celltype, "_TFreg_score_per_sample.RDS"))
# # TF在不同肿瘤中的鉴定情况
# # 要不要提前把全是0的TF删掉
# # 每个样本中的correlation怎么保留，要看P吗？
# # 之前的cor不看单个的P？
# 
# 
# TFreg <- lapply(unique(sam_info$Dataset), function(tumor) {
#     return(GSEABase::getGmt(paste0("res/02_SCENIC/regulon/", tumor, "_reg_metab.gmt")))
# })
# names(TFreg) <- unique(sam_info$Dataset)
# 
# 
# #### arrange path scores
# path_score <- lapply(unique(sam_info$Dataset), function(tumor) {
#     metab_pas <- data.frame(data.table::fread(paste0("res/02_SCENIC/aucell/", tumor, "_metab_aucell.csv")),
#                             row.names = 1, check.names = F)
#     used_sams = sam_info$samID[sam_info$Dataset == tumor]
#     tmp_res <- lapply(used_sams, function(sam) {
#         metab_pas[sam_cells[[sam]], ]
#     })
#     names(tmp_res) <- used_sams
#     return(tmp_res)
# })
# path_score_per_sample = do.call(c, path_score)
# saveRDS(path_score_per_sample, file = paste0("res/02_NMF/MPs/", celltype, "_path_score_per_sample.RDS"))


# calculate correlation ============================================================================
all(names(MP_scores_per_sample) == names(sig_score_per_sample))
source("src/03_NMF/00_tools_score_cor_heatmap.R")
source("src/03_NMF/00_tools_score_cor_pooled_heatmap.R")


## MP x MP
score_cor_pooled(score1 = lapply(MP_scores_per_sample, function(x) x+matrix(rnorm(nrow(x)*ncol(x), sd = 10e-5), nrow = nrow(x))), 
                  score2 = MP_scores_per_sample, sam_info = sam_info, max_cor = 1,
                  out2pdf = TRUE, fout = paste0("plot/02_NMF/", celltype, "_MPs_x_MPs_heatmap.pdf"))

## MPs x Gavish
score_cor_pooled(score1 = MP_scores_per_sample, score2 = Gavish_scores_per_sample,
                  sam_info = sam_info, max_cor = 1, 
                  out2pdf = TRUE, fout = paste0("plot/02_NMF/", celltype, "_Gavish_x_MPs_heatmap.pdf"), 
                  width = 18, height = 20)

## MP x Sig
score_cor_heatmap(score1 = MP_scores_per_sample, score2 = sig_score_per_sample, sam_info = sam_info, min_datasets = 4,
                  out2pdf = FALSE, fout = paste0("plot/02_NMF/", celltype, "_MPs_x_Sigs_heatmap.pdf"), max_cor = 1)

Imp_sig <- c("CD8.c01(Tn)", "CD8.c02(IL7R+ Tm)", "CD8.c05(GZMK+ early Tem)",
             "CD8.c06(GZMK+ Tem)", "CD8.c10(ZNF683+CXCR6+ Trm)", "CD8.c12(terminal Tex)")
Imp_sig <- c("CD4.c01(Tn)", "CD4.c13(Temra)", "CD4.c17(IFNG+ Tfh/Th1)", "CD4.c20(TNFRSF9+ Treg)")
score_cor_heatmap(score1 = MP_scores_per_sample, score2 = lapply(sig_score_per_sample, function(x) x[, Imp_sig]),
                  sam_info = sam_info, filter_path = TRUE, min_datasets = 4,
                  fout = paste0("plot/2.5_NMF/", celltype, "_MPs_x_Imp_Sigs_heatmap.pdf"))

## Path x Sig
score_cor_heatmap(score1 = path_score_per_sample, score2 = sig_score_per_sample, 
                  sam_info = sam_info, filter_path = TRUE, min_datasets = 4, 
                  fout = paste0("plot/2.5_NMF/", celltype, "_Pathways_x_Sigs_heatmap.pdf"))

Imp_sig <- c("CD8.c01(Tn)", "CD8.c02(IL7R+ Tm)", "CD8.c05(GZMK+ early Tem)",
             "CD8.c06(GZMK+ Tem)", "CD8.c10(ZNF683+CXCR6+ Trm)", "CD8.c12(terminal Tex)")
Imp_sig <- c("CD4.c01(Tn)", "CD4.c13(Temra)", "CD4.c17(IFNG+ Tfh/Th1)", "CD4.c20(TNFRSF9+ Treg)")
score_cor_heatmap(score1 = path_score_per_sample, score2 = lapply(sig_score_per_sample, function(x) x[, Imp_sig]), 
                  sam_info = sam_info, filter_path = TRUE, min_datasets = 4, 
                  fout = paste0("plot/2.5_NMF/", celltype, "_Pathways_x_Imp_Sigs_heatmap.pdf"))

## MP x Path
score_cor_heatmap(score1 = path_score_per_sample, score2 = MP_scores_per_sample, 
                  sam_info = sam_info, filter_path = TRUE, min_datasets = 4, 
                  out2pdf = TRUE, fout = paste0("plot/02_NMF/", celltype, "_Pathways_x_MPs_heatmap.pdf"), 
                  width = 18, height = 20)


## MP x TF
source("src/03_NMF/00_tools_score_cor_heatmap_TF.R")
TF_data = score_cor_heatmap_TF(score1 = TFreg_score_per_sample, score2 = MP_scores_per_sample, sam_info = sam_info, 
                     filter_path = TRUE, min_datasets = 4, width = 18, height = 20, max_cor = 1, 
                     fout = paste0("plot/02_NMF/", celltype, "/", celltype, "_TFreg_x_MPs_heatmap.pdf"))
write.table(TF_data, file = "res/02_NMF/TF_of_MMPs.tsv", quote = F, sep = "\t")



# II. CD8 T cells sig x sig
rm(list=ls())
source("src/02_NMF/00_tools_score_cor_pooled_heatmap.R")
sig_score_per_sample <- readRDS(file = paste0("res/02_NMF/MPs/CD8_T_sig_score_per_sample.RDS"))
sam_info_list <- readRDS("data/rds/NMF/cells_info_per_sam.RDS")$CD8_T
attach(sam_info_list)
MP_list <- readRDS("res/02_NMF/MPs/MP_list.RDS")$CD8_T

score_cor_pooled(score1 = lapply(sig_score_per_sample, function(x) x+matrix(rnorm(nrow(x)*ncol(x), sd = 10e-5), nrow = nrow(x))), 
                 score2 = sig_score_per_sample,
                 sam_info = sam_info, max_cor = 1, 
                 out2pdf = TRUE, fout = paste0("plot/02_NMF/", celltype, "_sig_x_sig_heatmap.pdf"), 
                 width = 18, height = 20)
sig_score_per_sample <- lapply(sig_score_per_sample, function(x) x[, c(1, 2, 5, 6, 10, 12)])
score_cor_pooled(score1 = lapply(sig_score_per_sample, function(x) x+matrix(rnorm(nrow(x)*ncol(x), sd = 10e-5), nrow = nrow(x))), 
                 score2 = sig_score_per_sample,
                 cluster_rows = FALSE, cluster_columns = FALSE, 
                 sam_info = sam_info, max_cor = 1, 
                 out2pdf = TRUE, fout = paste0("plot/02_NMF/", celltype, "_sig_x_sig_exhaustion_heatmap.pdf"), 
                 width = 18, height = 20)

CD8T_expr <- readRDS("data/rds/NMF/log_transformed/CD8_T_log2_all_samples.RDS")
oxphos_gene <- c(MP_list$`MMP3 OXPHOS1`, MP_list$`MMP8 OXPHOS2`)
CD8T_oxphos <- lapply(CD8T_expr, function(x) t(as.matrix(x[oxphos_gene, ])))
MP_scores_per_sample <- readRDS(paste0("res/02_NMF/MPs/CD8_T_MP_scores_p30_.7_.2_final_MinScore_0.8.RDS"))


score_cor_pooled(score1 = CD8T_oxphos,  
                 score2 = sig_score_per_sample,
                 cluster_rows = FALSE, cluster_columns = FALSE, 
                 sam_info = sam_info, R_cutoff = 0, max_cor = 0.5, 
                 out2pdf = FALSE, fout = paste0("plot/02_NMF/", celltype, "_sig_x_sig_exhaustion_heatmap.pdf"), 
                 width = 18, height = 20)

score_cor_pooled(score1 = MP_scores_per_sample,  
                 score2 = sig_score_per_sample,
                 cluster_rows = FALSE, cluster_columns = FALSE, 
                 sam_info = sam_info, R_cutoff = 0, max_cor = 0.5, 
                 out2pdf = FALSE, fout = paste0("plot/02_NMF/", celltype, "_sig_x_sig_exhaustion_heatmap.pdf"), 
                 width = 18, height = 20)

                      