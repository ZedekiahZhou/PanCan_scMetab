# ==================================================================================================
# Author: Zhou Zhe
# Program: Calculate MP abundance (Plot context specifiity)
# Version: 1.0
# Date: Aug 11, 2023
# ==================================================================================================

rm(list=ls())
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)

### Define the following:
para = "p30_.7_.2_final"
celltype = "Mast_cells"

for (celltype in c("Malignant", "Myeloid", "B_cells", "CD4_T", "CD8_T", "Endothelial", "Fibroblasts")) {
    load(paste0(paste0("res/02_NMF/MPs/", celltype, "_MPs_p30_.7_.2_final.RData")))
    MP_list <- readRDS("res/02_NMF/MPs/MP_list.RDS")
    
    programs_info <- data.frame(names = colnames(nmf_programs_original))
    programs_info$tumor_type = sub("_.+", "", programs_info$names)
    programs_info$dataset = sapply(strsplit(programs_info$names, "_"), function(y) {
        if (length(grep("valid", y)) > 0) paste(y[1:2], collapse = "_") else y[1]
    })
    cancer_list = split(programs_info$names, programs_info$tumor_type)
    
    abundance = plyr::ldply(1:length(Cluster_list), function(i) {
        MMP = names(MP_list[[celltype]])[i]
        plyr::ldply(names(cancer_list), function(cancer) {
            x = Cluster_list[[i]]  ### MP-related programs
            y = cancer_list[[cancer]]   ### cancer related programs
            observed = length(intersect(x, y))  ### observed number of MP-related programs in a cancer type
            expected = length(x)*length(y)/length(unlist(cancer_list))  ### expected number of MP-related programs in a cancer type
            A = log2((observed+1)/(expected+1))
            p = phyper(observed-1, length(y), 
                       length(unlist(cancer_list))-length(y), length(x), lower.tail = F)
            return(data.frame(Cluster = MMP, Cancer = cancer, 
                              observed = observed, expected = expected, 
                              MP_size = length(x), cancer_size = length(y), 
                              A = A, p_value = p))
        })
    })
    
    abundance$p_adj = p.adjust(abundance$p_value, method = "bonferroni")
    abundance = mutate(abundance, 
                       classification = ifelse((observed > 10) | (A > 1), 
                                               ifelse(p_adj < 0.05, "High_significant", "High"), 
                                               ifelse(((observed <= 10) & (observed >= 2)) | ((A <= 1) & (A > 0)), "Medium", 
                                                      ifelse((observed == 1) & ((A <= 0) & (A > -1.5)), "Low", 
                                                             ifelse(observed == 0, "Absent", "Absent"))))) 
    abundance$classification <- factor(abundance$classification, levels = c("Absent", "Low", "Medium", "High", "High_significant"))
    tmp <- abundance %>% group_by(Cluster) %>% summarise(n = sum(classification != "Absent")) %>% arrange(desc(n))
    abundance$Cluster <- factor(abundance$Cluster, levels = tmp$Cluster)
    
    my_cols = c("white", viridis::magma(4, direction = -1, begin = 0.18, end = 1))
    names(my_cols) <- c("Absent", "Low", "Medium", "High", "High_significant")
    
    # tmp.text <- element_text(family="sans", size=8)
    # ggplot(data = abundance, aes(x=Cancer, y=Cluster, fill=classification, color=classification)) + 
    #     geom_tile() + 
    #     scale_color_manual(values = my_cols) + 
    #     scale_fill_manual(values = my_cols) + 
    #     theme_bw() +
    #     scale_x_discrete(expand = c(0,0)) + 
    #     scale_y_discrete(expand = c(0,0)) + 
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
    #           axis.ticks = element_blank(), text = tmp.text, 
    #           legend.key = element_rect(colour = "black", linewidth = 0.5), 
    #           legend.key.size = unit(4, "mm"))
    # ggsave(paste0("plot/2.5_NMF/", celltype, "_MPs_context_specificity_tumor_type_p30_.7_.2_final.pdf"), 
    #        width = 4, height = 2)
    
    plot_data = reshape2::dcast(abundance, Cluster ~ Cancer, value.var = "classification")
    rownames(plot_data) = plot_data$Cluster
    tmp.text <- gpar(fontfamily="sans", fontsize = 8)
    print(celltype)
    pdf(paste0("plot/02_NMF/", celltype, "_MPs_context_specificity_tumor_type.pdf"), 
        width = 7, height = 7)
    print(Heatmap(plot_data[rev(1:nrow(plot_data)), -1], col = my_cols, border = T, 
            row_names_gp = tmp.text, column_names_gp = tmp.text, row_names_side = "left", 
            row_title_gp = tmp.text, column_title_gp = tmp.text, 
            height = nrow(plot_data) * unit(0.3, "cm"), 
            width = ncol(plot_data) * unit(0.3, "cm"), 
            heatmap_legend_param = list(border = "black"), 
            name = " ", 
            cluster_rows = F, cluster_columns = F))
    dev.off()
    
}
