rm(list = ls())
library(Seurat)
library(scMetab)
library(RColorBrewer)
library(ggplot2)

seed.use = 2021
tmp.text <- element_text(family="sans", size=18)

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

seu.list <- lapply(1:10, function(i) {
    tumor <- dataset$DataSets[i]
    print(tumor)
    seu_obj <- readRDS(paste0("data/rds/PAS/", tumor, ".TpN.PAS.unintegrated.rds"))
    
    set.seed(seed.use)
    seu_obj <- seu_obj[, sample(x = colnames(seu_obj), size = 0.1*dim(seu_obj)[2])]
    seu_obj$dataset <- tumor
    
    gc()
    return(seu_obj)
})
set.seed(NULL)
names(seu.list) <- dataset$DataSets
seu.list <- seu.list[sort(names(seu.list))]

seu <- merge(seu.list[[1]], seu.list[-1])
seu$celltype_x_TorN

## plot malignant cells 
malig <- subset(seu, celltype_x_TorN %in% c("Malignant_T", "Epithelial_N") & dataset != "BRCA_valid1")

tmp.text <- element_text(family="sans", size=24)
pdf("plot/2_TvN/DE/Malignant_VlnPlot.pdf", width = 24, height = 9)
DefaultAssay(malig) <- "RNA"
VlnPlot(malig, features = c("CA2", "GAPDH", "LDHA", "ATP1A1", "SCD", "ASAH1"), pt.size = 0, 
        group.by = "dataset", split.by = "TorN", cols = TvN_color, ncol = 4) & 
    geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
    ylab("") & xlab("") &
    theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, axis.title = tmp.text)
DefaultAssay(malig) <- "Metab"
VlnPlot(malig, features = c("Fatty Acid Biosynthesis", "Citric Acid Cycle", 
                            "Oxidative Phosphorylation", "Glycolysis and Gluconeogenesis", 
                            "Histidine Metabolism", "Creatine Metabolism"), 
        pt.size = 0, group.by = "dataset", split.by = "TorN", cols = TvN_color, ncol = 4) & 
    geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
    ylab("") & xlab("") &
    theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, axis.title = tmp.text)
dev.off()


## plot myeloid cells 
myeloid <- subset(seu, celltype == "Myeloid" & dataset != "BRCA_valid1")
pdf("plot/2_TvN/DE/Myeloid_VlnPlot.pdf", width = 24, height = 4.5)
DefaultAssay(myeloid) <- "RNA"
VlnPlot(myeloid, features = c("FTH1", "SLC6A6", "GAPDH"), pt.size = 0, 
        group.by = "dataset", split.by = "TorN", cols = TvN_color, ncol = 4) & 
    geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
    ylab("") & xlab("") &
    theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, axis.title = tmp.text)
DefaultAssay(myeloid) <- "Metab"
VlnPlot(myeloid, features = c("Oxidative Phosphorylation", "Glycolysis and Gluconeogenesis", 
                            "Histidine Metabolism", "Lysine Metabolism"), 
        pt.size = 0, group.by = "dataset", split.by = "TorN", cols = TvN_color, ncol = 4) & 
    geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) &
    ylab("") & xlab("") &
    theme(text=tmp.text, legend.text = tmp.text, axis.text = tmp.text, axis.title = tmp.text)
dev.off()





