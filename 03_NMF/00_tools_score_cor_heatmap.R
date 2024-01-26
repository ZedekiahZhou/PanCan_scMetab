# -------------------------------------------------------------------------------------------
# Function for heatmap between different scores: 
#   correlations calculated within each tumor, and then averaged across tumors from the same tumor type
# ------------------------------------------------------------------------------------------- 

score_cor_heatmap <- function(score1, # a list contains scores for cells of each sample
                              score2, 
                              sam_info, # a data.frame contains samID, Tumor_Type and Dataset
                              p_cutoff = 0.05, R_cutoff = 0.1, 
                              filter_path = FALSE, # wheather to filter path (rows)
                              min_datasets = 3, # minimal consistent datasets (in at least a Sig_path) to keep a pathway 
                              out2pdf = TRUE, # wheather to output to *.pdf file
                              fout, # output file name (must be *.pdf)
                              width = 18, height = 9,
                              max_cor = 0.5
) {
    require(ComplexHeatmap)
    pearson_cor <- mapply(function(x, y) {
        tmp_res <- reshape2::melt(cor(x, y, method = "pearson"))
        colnames(tmp_res) <- c("MMPs", "Sigs", "R")
        tmp_res$ID <- paste(tmp_res$MMPs, tmp_res$Sigs, sep = "_")
        return(tmp_res)
    }, score1, score2, SIMPLIFY = F)
    
    spearman_cor <- mapply(function(x, y) {
        tmp_res <- reshape2::melt(cor(x, y, method = "spearman"))
        colnames(tmp_res) <- c("MMPs", "Sigs", "Rho")
        tmp_res$ID <- paste(tmp_res$MMPs, tmp_res$Sigs, sep = "_")
        return(tmp_res)
    }, score1, score2, SIMPLIFY = F)
    
    cor_df_list <- lapply(unique(sam_info$Dataset), function(tumor) {
        sams = sam_info$samID[sam_info$Dataset == tumor]
        cor_df = plyr::adply(pearson_cor[[sams[1]]][c("ID", "MMPs", "Sigs")], 1, function(x) {
            all_pearson <- sapply(pearson_cor[sams], function(y) y$R[y$ID == x$ID])
            all_spearman <- sapply(spearman_cor[sams], function(y) y$Rho[y$ID == x$ID])
            data.frame(ID = x$ID, MMPs = as.character(x$MMPs), Sigs = as.character(x$Sigs), 
                       Mean_R = mean(all_pearson, na.rm = T), P_pearson = t.test(all_pearson)$p.value, 
                       Mean_Rho = mean(all_spearman, na.rm = T), P_spearman = t.test(all_spearman)$p.value)
        })
        cor_df$Padj_pearson = p.adjust(cor_df$P_pearson, method = "BH")
        cor_df$Padj_spearman = p.adjust(cor_df$P_spearman, method = "BH")
        cor_df$R = ifelse(cor_df$Padj_pearson >= p_cutoff | cor_df$Padj_spearman >= p_cutoff | abs(cor_df$Mean_R) < R_cutoff, 
                          0, cor_df$Mean_R)
        cor_mtx <- reshape2::acast(cor_df, MMPs ~ Sigs, fun.aggregate = mean, value.var = "R")
        colnames(cor_mtx) <- paste(colnames(cor_mtx), tumor, sep = "_")
        return(cor_mtx)
    })
    names(cor_df_list) <- unique(sam_info$Dataset)
    cor_df_comnbined <- data.frame(do.call(cbind, cor_df_list), check.names = F)
    
    Sig_path = sort(colnames(score2[[1]]))
    column_annot <- rep(unique(sam_info$Dataset), each = length(Sig_path))
    column_split <- rep(Sig_path, times = length(unique(sam_info$Dataset)))
    
    if (filter_path) {
        n_significance <- sapply(split.default(cor_df_comnbined, column_split), function(x) rowSums((!is.na(x)) & (x != 0)))
        idx <- rowSums(n_significance >= min_datasets) > 0
        plot_data = cor_df_comnbined[idx, ]
    } else {
        plot_data = cor_df_comnbined
    }
    
    
    tmp.text <- gpar(fontfamily="sans", fontsize = 8)
    if (max_cor == 0.5) {
        myCols <- rev(brewer.pal(n=5, name="RdBu"))
    } else {
        myCols <- rev(brewer.pal(n=11, name="RdBu"))[c(1, 2, 6, 10, 11)]
    }
    plot_data[is.na(plot_data)] = 0
    
    if (out2pdf) pdf(fout, width = width, height = height)
    p = Heatmap(plot_data, col = circlize::colorRamp2(c(-max_cor, -max_cor/2, 0, max_cor/2, max_cor), myCols), 
            rect_gp = gpar(col = "white"), 
            cluster_columns = FALSE, cluster_rows = TRUE, clustering_method_rows = "ward.D2", 
            row_names_side = "left", show_row_dend = FALSE, 
            row_names_max_width = unit(10, "cm"),
            row_names_gp = tmp.text, column_names_gp = tmp.text, 
            row_title_gp = tmp.text, column_title_gp = tmp.text, 
            column_title_rot = 90, 
            height = unit(0.282*nrow(plot_data), "cm"), 
            width = unit(0.15*ncol(plot_data), "cm"), 
            heatmap_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
            top_annotation = HeatmapAnnotation(Tumor = column_annot, col = list(Tumor = dataset_color), 
                                               simple_anno_size = unit(2, "mm"), 
                                               show_annotation_name = F, 
                                               annotation_legend_param = list(title_gp = tmp.text, labels_gp = tmp.text), 
                                               which = "column"), 
            column_split = column_split, gap = unit(1, "mm"), 
            show_column_names = FALSE, column_names_side = "top", 
            column_names_max_height = unit(6, "cm"),
            border = TRUE, name = "Correlation")
    print(p)
    if (out2pdf) dev.off()
}
