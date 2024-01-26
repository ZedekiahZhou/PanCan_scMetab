# =================================================================================================#
# Author: Zhou Zhe
# Program: 
#   1. Annotate MPs manually
#   2. Heatmap of MP
#   3. plot MP scores for all MP genes across NMF programs
#   4. The realationship between MMP3 and MMP12
# 
# Version: 1.0
# Date: Aug 11, 2023
# =================================================================================================#

rm(list=ls())
library(ComplexHeatmap)
all_celltypes <- c("Malignant", "Myeloid", "B_cells", "CD4_T", "CD8_T", "Mast_cells", "Endothelial", "Fibroblasts")
ct = "Malignant"

# I. Annotate =====================================
## Malignant cells
load(file = paste0("res/02_NMF/MPs/Malignant_MPs_p30_.7_.2_final.RData"))
anno <- c("MMP1 DNA Synthesis", "MMP2 Steroid Metabolism", "MMP3 OXPHOS1", 
          "MMP4 Glycolysis", "MMP5 CAC", "MMP6 Transport", "MMP7", "MMP8 Glycosylation", 
          "MMP9", "MMP10 AA Metabolism", "MMP11", "MMP12 OXPHOS2", "MMP13", 
          "MMP14", "MMP15")
Malignant_MP <- MP_list
names(Malignant_MP) <- anno

## Myeloid 
load(file = paste0("res/02_NMF/MPs/Myeloid_MPs_p30_.7_.2_final.RData"))
anno <- c("MMP1", "MMP2", "MMP3 DNA Synthesis", 
          "MMP4", "MMP5 OXPHOS", "MMP6 Eicosanoid Metabolism", 
          "MMP7 CAC", "MMP8", "MMP9 Glycolysis", "MMP10", "MMP11 Fatty Acid Biosynthesis", 
          "MMP12", "MMP13 Transport", "MMP14", "MMP15", "MMP16", "MMP17")
Myeloid_MP <- MP_list
names(Myeloid_MP) <- anno
a = sapply(Myeloid_MP, function(x) sapply(Malignant_MP, function(y) length(intersect(x, y))))

## B_cells
load(file = paste0("res/02_NMF/MPs/B_cells_MPs_p30_.7_.2_final.RData"))
anno <- c("MMP1 Glycosylation", "MMP2", "MMP3 DNA Synthesis", 
          "MMP4 AA Metabolism", "MMP5", "MMP6 OXPHOS", 
          "MMP7", "MMP8", "MMP9", "MMP10")
B_cells_MP <- MP_list
names(B_cells_MP) <- anno
a = sapply(B_cells_MP, function(x) sapply(Malignant_MP, function(y) length(intersect(x, y))))

## CD4_T
load(file = paste0("res/02_NMF/MPs/CD4_T_MPs_p30_.7_.2_final.RData"))
anno <- c("MMP1 DNA Synthesis", "MMP2 Transport", "MMP3", 
          "MMP4 OXPHOS", "MMP5", "MMP6")
CD4_T_MP <- MP_list
names(CD4_T_MP) <- anno
a = sapply(CD4_T_MP, function(x) sapply(Myeloid_MP, function(y) length(intersect(x, y))))

## CD8_T
load(file = paste0("res/02_NMF/MPs/CD8_T_MPs_p30_.7_.2_final.RData"))
anno <- anno <- c("MMP1 DNA Synthesis", "MMP2 Transport", "MMP3 OXPHOS1", 
                  "MMP4", "MMP5", "MMP6", 
                  "MMP7", "MMP8 OXPHOS2", "MMP9", "MMP10", "MMP11", "MMP12", "MMP13")
CD8_T_MP <- MP_list
names(CD8_T_MP) <- anno
a = sapply(CD8_T_MP, function(x) sapply(CD4_T_MP, function(y) length(intersect(x, y))))

## Endothelial
load(file = paste0("res/02_NMF/MPs/Endothelial_MPs_p30_.7_.2_final.RData"))
anno <- c("MMP1", "MMP2", "MMP3", "MMP4", "MMP5 OXPHOS",
          "MMP6 CAC", "MMP7", "MMP8", "MMP9", "MMP10", "MMP11", "MMP12")
Endothelial_MP <- MP_list
names(Endothelial_MP) <- anno
a = sapply(Endothelial_MP, function(x) sapply(Malignant_MP, function(y) length(intersect(x, y))))

## Fibroblasts
load(file = paste0("res/02_NMF/MPs/Fibroblasts_MPs_p30_.7_.2_final.RData"))
anno <- c("MMP1", "MMP2", "MMP3", 
          "MMP4", "MMP5", "MMP6 DNA Synthesis", 
          "MMP7 OXPHOS", "MMP8", "MMP9", "MMP10 Glycolysis", 
          "MMP11", "MMP12", "MMP13", "MMP14", "MMP15", "MMP16")
Fibroblasts_MP <- MP_list
names(Fibroblasts_MP) <- anno
a = sapply(Fibroblasts_MP, function(x) sapply(Malignant_MP, function(y) length(intersect(x, y))))

MP_all <- list(Malignant = Malignant_MP, Myeloid = Myeloid_MP, B_cells = B_cells_MP, 
               CD4_T = CD4_T_MP, CD8_T = CD8_T_MP, Endothelial = Endothelial_MP, Fibroblasts = Fibroblasts_MP) 
saveRDS(MP_all, file = "res/02_NMF/MPs/MP_list.RDS")

writexl::write_xlsx(MP_all, path = "res/02_NMF/MPs/MP_list.xlsx")

## Jaccard index between Malignant MPs and TME MPs
tmp.text <- gpar(fontfamily="sans", fontsize = 8)
pdf("plot/02_NMF/Jaccard_similarity_Malignant_x_TME.pdf", width = 10, height = 10)
jaccard <- lapply(names(MP_all)[-1], function(MP) {
    res = sapply(MP_all$Malignant, function(x) {
        sapply(MP_all[[MP]], function(y) length(intersect(x, y))/length(union(x, y)))
    })
    print(Heatmap(res, col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                  show_row_dend = F, show_column_dend = F, 
                  height = nrow(res) * unit(0.3, "cm"), 
                  width = ncol(res) * unit(0.3, "cm"), 
                  row_title = paste0(MP, "_MPs"), column_title = "Malignant_MPs", 
                  row_title_gp = tmp.text, column_title_gp = tmp.text, 
                  row_names_gp = tmp.text, column_names_gp = tmp.text, 
                  #rect_gp = gpar(col = "white"), gap = unit(1, "mm"), 
                  name = MP))
})
names(jaccard) = names(MP_all)[-1]
dev.off()


# II. replot Malignant cells include only MPs ================================
####  Sort Jaccard similarity plot according to new clusters:
load(file = paste0("res/02_NMF/MPs/Malignant_MPs_p30_.7_.2_final.RData"))
source("src/02_NMF/custom_magma.R")
inds_sorted <- c()

for (j in 1:length(Cluster_list)){
    
    inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))
    
}
#inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered
inds_new <- inds_sorted ## remove unclustered programs

nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new]) 
nmf_intersect_meltI_NEW$Var2 <- factor(nmf_intersect_meltI_NEW$Var2, levels = rev(levels(nmf_intersect_meltI_NEW$Var1)))

ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
    geom_raster() + 
    scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
    scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
    theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
ggsave(filename = paste0("plot/02_NMF/", ct, "_MP_cluster_p30_7_.2_Clustered_only.pdf"), width = 10, height = 9)


# III. NMF scores plot (not used) =====================
load("res/02_NMF/MPs/Malignant_MPs_p30_.7_.2_final.RData")
Genes_nmf_w_basis <- readRDS("res/02_NMF/NMF/malig/Genes_nmf_w_basis.RDS")
Genes_nmf_w_basis_robust <- lapply(Genes_nmf_w_basis, function(x) {
    y = data.frame(x[, colnames(x) %in% colnames(nmf_programs_original)], check.names = F)
    y$Gene = rownames(y)
    return(y)
})
Genes_nmf_w_basis_robust = Reduce(function(x, y) merge(x, y, by = "Gene", all = T), 
                                  Genes_nmf_w_basis_robust)
rownames(Genes_nmf_w_basis_robust) = Genes_nmf_w_basis_robust$Gene
Genes_nmf_w_basis_robust <- Genes_nmf_w_basis_robust[-1]
Genes_nmf_w_basis_robust <- Genes_nmf_w_basis_robust[rev(unlist(MP_list)), unlist(Cluster_list)]
Genes_nmf_w_basis_robust <- data.frame(t(Genes_nmf_w_basis_robust), check.names = F)

plot_data = reshape2::melt(as.matrix(Genes_nmf_w_basis_robust))
my_color = brewer.pal(9, "GnBu")[c(1,2,5,6,8,9)]
br = c(0, 0.25, 0.26, 0.49, 0.5, 1)
br = c(0, 0.49, 0.5, 0.75, 0.76, 1)
ggplot(plot_data, aes(x=Var1, y=Var2, fill=value, color=value)) + 
    geom_raster() + 
    scale_color_gradientn(limits=c(0, 0.2), colours = my_color, 
                          values = br, na.value = "white", guide = "legend") + 
    scale_fill_gradientn(limits=c(0, 0.2), colours = my_color, 
                         values = br, na.value = "white", guide = "colourbar") + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

# a = data.frame(sapply(split(Genes_nmf_w_basis_robust, rep(1:15, each = 30)), function(x) rowMeans(x, na.rm = T)))
# a$MP = apply(a, 1, which.max)

# Genes_nmf_w_basis_robust <- do.call(cbind, Genes_nmf_w_basis_robust)
MP_genes <- unlist(MP_list)
n_list <- lapply(Genes_nmf_w_basis_robust, function(x) intersect(rownames(x), MP_genes))
