# ==================================================================================================
# Author: Zhou Zhe
# Program: Program to plot Malignant MPs gene heatmap for all cells in each sample 
#          (Modified from Gavish et. al.)
# Version: 1.0
# Date: Aug 11, 2023
# ==================================================================================================

# ------------------------------------------------------------------------------------------- 

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(scales)
library(Matrix)
library(parallel)

### Define the following:
celltype = "Malignant"
My_study_sparse  = readRDS(paste0("data/rds/NMF/log_transformed/", sub(" ", "_", celltype), "_log2_all_samples.RDS"))
used_MP = "MMP_7"

para = "p30_.7_.2_final"
MinScore = 0.8
load(paste0("res/02_NMF/MPs/", celltype, "_MPs_", para, ".RData"))
MP_all <- readRDS("res/02_NMF/MPs/MP_list.RDS")
heatCols  = readRDS("src/02_NMF/heatCols.RDS.gz")  # can be found in the Github repository
MP_score <- readRDS("res/02_NMF/MPs/Malignant_MP_scores_p30_.7_.2_final_MinScore_0.8.RDS")

### plot function
color.scheme <- colorRampPalette(c(heatCols))(n=333)
Plot_heatmap <- function(sam) { 
    M <- L_plot2[[sam]]
    M_new2        <- apply(M, 2, rev)
    M_meltII      <-  reshape2::melt(t(M_new2)) 
    M_meltII$Var2 <- factor(M_meltII$Var2)
    
    G <- ggplot(data = M_meltII, aes(x=Var1, y=Var2, fill=value, color=value)) + 
        geom_raster() + 
        scale_color_gradient2(limits=c(-4,4), low=color.scheme[1:111],  mid =color.scheme[112:222], high = color.scheme[223:333], midpoint = 0, oob=squish, name=NULL) +                                
        scale_fill_gradient2(limits=c(-4,4), low=color.scheme[1:111],  mid =color.scheme[112:222], high = color.scheme[223:333], midpoint = 0, oob=squish, name=NULL)  +
        scale_x_discrete(labels = rep("", ncol(M))) + 
        theme(panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 8), legend.title = element_text(size=8), legend.text = element_text(size = 8), legend.text.align = 0.5, legend.justification = "bottom" ) + 
        theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
        theme(axis.title.y=element_blank(), axis.ticks.y=element_blank()) +
        ggtitle(sam)
    
    return(G)
}

Plot_box <- function(sam) {
    M <- reshape2::melt(t(apply(L_plot[[sam]], 2, rev)))
    ggplot(M, aes(x = Var2, y = value)) + 
        geom_boxplot(outlier.shape = NA) + 
        coord_flip() + 
        theme_classic() + 
        theme(axis.title = element_blank(), axis.text.y=element_blank()) + 
        ggtitle(" ")
                        
}



for (used_MP in names(MP_list)) {
    ### select MP genes matrix and centering
    MP_genes <- MP_list[[used_MP]]
    L_plot        <- lapply(seq_along(My_study_sparse), function(I) {
        idx = match(MP_genes , rownames(My_study_sparse[[I]]))
        res = as.matrix(My_study_sparse[[I]][idx[!is.na(idx)] , ])
    }) 
    names(L_plot) <- names(My_study_sparse)
    
    L_plot2 <- lapply(L_plot, function(x) {
        res = x - rowMeans(x)
        res = res[, order(colSums(res))]
    })
    
    P1 <- lapply(names(L_plot2), Plot_heatmap)
    P2 <- lapply(names(L_plot), Plot_box)
    # P_list <- do.call(c, lapply(1:length(P1), function(i) c(P1[i], P2[i])))
    # P_list <- do.call(ggarrange, c(P_list, ncol=1, nrow = 2))
    
    pdf(paste0("plot/02_NMF/", celltype, "_", used_MP, "_gene_heatmap_raster.pdf"), width = 15, height = 5)
    lapply(1:length(P1), function(i) {print(ggarrange(P1[[i]], P2[[i]], nrow = 1, widths = c(4, 1)))})
    dev.off() 
}
