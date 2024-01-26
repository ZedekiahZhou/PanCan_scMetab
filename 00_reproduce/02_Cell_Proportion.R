# ==================================================================================================
# Author: Zhou Zhe
# Program: Plot Sample Counts, Cell Counts & Celltype proportion
# Version: 1.0
# Date: Aug 12, 2022
# ==================================================================================================
library(ggplot2)
library(ggpubr)
library(scMetab)
library(tidyverse)

dir.used <- "00_reproduce/AllGene_clustering/"
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

# I. Cell type proportions ------------------------------------------------------------------------
celltype_counts.list <- lapply(dataset$DataSets, function(tumor) {
    res = read.csv(file = paste0("res/", dir.used, tumor, "_celltype_counts.csv"), 
                   colClasses = c("character", "character", "integer"))
    res$dataset <- tumor
    res
})
celltype_counts_pooled <- do.call(rbind, celltype_counts.list) 


tmpfile <- paste0("plot/", dir.used, "Cell_Proportion.pdf")
tmp.text <- element_text(family="sans", size=8)
pdf(file = tmpfile, width = 6, height = 3)
ggplot(data = celltype_counts_pooled, mapping = aes(x = Group, y = Freq, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    xlab("Group") + ylab("Proportion") +
    scale_fill_manual(values = celltype_color) + 
    theme_classic() + 
    theme(axis.line=element_line(size = 0.3, colour="black"),
          legend.key.size = unit(0.4, "cm"),
          strip.text = tmp.text, 
          text = tmp.text) + 
    facet_wrap(~ dataset, nrow = 1)
ggbarplot(celltype_counts_pooled, x = "Group", y = "Freq", fill = "CellType", size = 0.3, 
          legend = "right", position = position_fill()) +
    scale_fill_manual(values = celltype_color) +
    theme(axis.line=element_line(size = 0.3, colour="black"),
          legend.key.size = unit(0.4, "cm"), 
          strip.text = tmp.text, 
          text = tmp.text) +
    ylab("Proportion") +
    facet_wrap(~ dataset, nrow = 1)
dev.off()


tmpfile <- paste0("plot/", dir.used, "Cell_Proportion_vertical.pdf")
tmp.text <- element_text(family="sans", size=6)
pdf(file = tmpfile, width = 3.5, height = 4)
ggplot(data = celltype_counts_pooled, mapping = aes(x = Group, y = Freq, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    xlab("Group") + ylab("Proportion") +
    scale_fill_manual(values = celltype_color) + 
    theme_classic() + 
    theme(axis.line=element_line(size = 0.3, colour="black"),
          legend.key.size = unit(0.4, "cm"),
          strip.text = tmp.text, 
          text = tmp.text) + 
    coord_flip() + 
    facet_wrap(~ dataset, ncol = 1, strip.position = "right", dir = "v")
dev.off()



## II. Samples and Cells --------------------------------------------------------------------------
sam.list <- lapply(dataset$DataSets, function(tumor) {
    res = read.csv(file = paste0("res/", dir.used, tumor, "_nSam_nCell_of_TandN.csv"))
    if (tumor == "BRCA_valid1") res$N = 0
    res$dataset <- tumor
    res
})
samDf <- do.call(rbind, sam.list)

nSam <- samDf %>%
    group_by(dataset) %>%
    summarise(N = sum(N > 0), T = sum(T > 0)) %>%
    reshape2::melt(id = 1, measure.vars = 2:3)

tmp.text <- element_text(family="sans", size=8)
p1 <- ggbarplot(nSam, x = "dataset", y = "value", fill = "variable", size = 0.3) +
    scale_fill_manual(values = TvN_color) + 
    xlab("Dataset") + ylab("Number of Samples") + 
    theme(axis.line=element_line(size = 0.3, colour="black"),
          legend.key.size = unit(0.4, "cm"), 
          text = tmp.text, axis.text.x = element_text(angle = 45, hjust = 1))

nCells <- samDf %>%
    group_by(dataset) %>%
    summarise(N = sum(N), T = sum(T)) %>%
    reshape2::melt(id = 1, measure.vars = 2:3)

p2 <- ggbarplot(nCells, x = "dataset", y = "value", fill = "variable", size = 0.3) +
    scale_fill_manual(values = TvN_color) + 
    xlab("Dataset") + ylab("Number of Cells") + 
    theme(axis.line=element_line(size = 0.3, colour="black"),
          legend.key.size = unit(0.4, "cm"), 
          text = tmp.text, axis.text.x = element_text(angle = 45, hjust = 1))

tmpfile <- paste0("plot/", dir.used, "Number_of_Samples.pdf")
ggsave(tmpfile, p1 + p2, width = 6, height = 3)
