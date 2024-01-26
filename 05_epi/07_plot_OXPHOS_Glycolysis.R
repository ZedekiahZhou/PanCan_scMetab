# ==================================================================================================
# Author: Zhou Zhe
# Program: Celltype TvN
# Version: 1.0
# Date: Jul 13, 2022
# ==================================================================================================

# init -------
rm(list=ls())

library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggsci)
library(scMetab)
library(RColorBrewer)
library(patchwork)
library(data.table)
# detach("package:scMetab", unload = TRUE)

n_cores <- 8
seed.use = 2021

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))

seu.list <- lapply(1:10, function(i) {
    tumor <- dataset$DataSets[i]
    print(paste0(tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(file = paste0("data/rds/PAS/", tumor, ".TpN.PAS.unintegrated.rds"))
    DefaultAssay(seu_obj) <- "Metab"
    seu_obj <- subset(seu_obj, TorN == "T" & celltype == "Malignant")
    seu_obj$dataset <- tumor
    
    gc()
    return(seu_obj)
})
names(seu.list) <- dataset$DataSets
seu.list <- seu.list[sort(names(seu.list))]


## density_plot function
density_plot <- function(path, 
                         path2 = NULL, 
                         dims = 1, 
                         assay = "Metab"
) {
    plot.list = lapply(1:nrow(dataset), function(i) {
        tumor <- dataset$DataSets[i]
        metab <- t(as.matrix(GetAssayData(seu.list[[tumor]], slot = "data", assay = assay)[path, ]))
        metab <- data.frame(metab, check.names = F)
        
        if (!is.null(path2)) {
            sig_pas <- data.frame(fread(paste0("res/5_SCENIC/aucell/", tumor, "_celltypeSig_aucell.csv")),
                              row.names = 1, check.names = F)[colnames(seu.list[[tumor]]), ]
            sig_pas <- sig_pas[, path2]
            metab <- cbind(metab, sig_pas)
        }
        
        # add noise to avoid error of geom_density_2d_filled()
        noise <- matrix(runif(nrow(metab)*ncol(metab), min = 0, max = 1e-10), ncol = ncol(metab))
        metab <- metab + noise
        
        if (dims == 1) {
            res.list <- lapply(colnames(metab), function(x) {
                ggplot(metab, aes(.data[[x]])) + 
                    geom_density(color = "black", fill = "grey") + 
                    theme_bw() + 
                    theme(legend.position = "none", panel.grid = element_blank())
            })
        } else if (dims == 2) {
            if (is.null(path2)) {
                df_comb = t(combn(path, m = 2))
            } else {
                df_comb = expand.grid(path, path2)
            }
            
            res.list <- apply(df_comb, 1, function(x) {
                ggplot(metab, aes(.data[[x[1]]], .data[[x[2]]])) + 
                    geom_density_2d_filled(contour_var = "ndensity") + 
                    scale_fill_manual(values = c("white", brewer.pal(9, "YlOrRd"))) + 
                    theme_bw() + 
                    theme(legend.position = "none", panel.grid = element_blank())
            })
        }
        
        res = wrap_plots(res.list, ncol = 4) + plot_annotation(tumor)
        return(res)
    })
    return(plot.list)
}

## 2D density
path.imp = c("Oxidative Phosphorylation", "Glycolysis and Gluconeogenesis", 
         "Fatty Acid Biosynthesis", "Citric Acid Cycle", 
         "Glutathione Metabolism", "Nucleotide Metabolism", 
         "Pentose Phosphate Pathway", "Glutamate metabolism")
path.lipid = c("Bile Acid Biosynthesis", 
         "Eicosanoid Metabolism", "Fatty Acid Biosynthesis", 
         "Fatty Acid Metabolism", "Fatty Acids Oxidation, mitochondrial", 
         "Fatty Acids Oxidation, peroxisomal", "Glycerophospholipid Metabolism",
         "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
         "Ketone Bodies Metabolism", "Sphingolipid Metabolism", 
         "Steroid Metabolism", "Triacylglycerol Synthesis")
cancer.path <- c("EMT", "Inflammatory response", "Angiogenesis", "Apoptosis", 
                 "mTORC1 signaling", "DNA repair", "G2M checkpoint")

# df = t(combn(path1, m = 2))

# 2D density for Imp path
pdf("plot/3_epi/Pathway_2D_density_Imp.pdf", width = 16, height = 28)
print(density_plot(path.imp, dims = 2))
dev.off()

pdf("plot/3_epi/Pathway_2D_density_Imp_x_hallmark.pdf", width = 16, height = 56)
print(density_plot(path = path.imp, path2 = cancer.path, dims = 2))
dev.off()

# 2D density for lipid path
pdf("plot/3_epi/Pathway_2D_density_lipid.pdf", width = 16, height = 78)
print(density_plot(path.lipid, dims = 2))
dev.off()

# 1D density for Imp path                                                                    
pdf("plot/3_epi/Pathway_1d_density.pdf", width = 16, height = 44)
print(density_plot(path = NULL, dims = 1))
dev.off()

# values = colorRampPalette(rev(brewer.pal(11, "RdGy")[1:6]))(10)



