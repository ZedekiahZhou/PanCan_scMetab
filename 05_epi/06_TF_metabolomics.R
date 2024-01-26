# ====================================================================
# Author: Zhou Zhe
# Function: Analysi metabolomics data of TF knockdown
# Version: 1.0
# Date: Mar 24, 2023
# ====================================================================

# 0. INIT 
rm(list = ls())
library(dplyr)
library(RColorBrewer)
library(ggplot2)
tmp.text <- element_text(family="sans", size=12)

# I. Prepare data =================
# sam_cols <- c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5), paste0("QC", 1:3))
# sam_groups <- c(rep(c("scramble", "HDAC2_kd", "FOSL1_kd"), each = 5), rep("QC", 3))

## Load H1299_HDAC2_p_FOSL1
input_dir = "data/metabolomics/H1299_HDAC2_p_FOSL1/"
header <- read.delim(paste0(input_dir, "header.txt"))
sam_cols <- c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5), paste0("QC", 1:3))

pos <- read.delim(paste0(input_dir, "Area_0_2022927942.txt"), skip = 4, check.names = F, quote = "", comment.char = "") %>%
    mutate(ID = paste0(.data$`Average Mz`, "__", .data$`Average Rt(min)`*60), 
           ID2 = paste0("POS", .data$`Alignment ID`), 
           Name = sub(";.+", "", sub('"$', '', sub('^"', '', .data$`Metabolite name`))))
idx <- intersect(grep("Unknown", pos$Name, invert = T), grep("^w/o", pos$Name, invert = T))
pos <- pos[idx, c("ID2", "Name", sam_cols)]

neg <- read.delim(paste0(input_dir, "Area_0_2022927951.txt"), skip = 4, check.names = F, quote = "", comment.char = "") %>%
    mutate(ID = paste0(.data$`Average Mz`, "__", .data$`Average Rt(min)`*60), 
           ID2 = paste0("NEG", .data$`Alignment ID`), 
           Name = sub(";.+", "", sub('"$', '', sub('^"', '', .data$`Metabolite name`))))  
idx <- intersect(grep("Unknown", neg$Name, invert = T), grep("^w/o", neg$Name, invert = T))
neg <- neg[idx, c("ID2", "Name", sam_cols)]

hdac2 <- rbind(pos, neg) 
hdac2$CV <- matrixStats::rowSds(as.matrix(hdac2[paste0("QC", 1:3)]))/
    rowMeans(as.matrix(hdac2[paste0("QC", 1:3)]))
hdac2 <- hdac2[!is.na(hdac2$CV) & hdac2$CV < 0.3, c(1:12, 21)]
hdac2_log <- log2(hdac2[, 3:12] + 1)
hdac2$log2FC <- log2(rowMeans(hdac2[paste0("B", 1:5)])/rowMeans(hdac2[paste0("A", 1:5)]))
tmp <- plyr::adply(hdac2_log, 1, function(x) {
    res = t.test(x[6:10], x[1:5])
    return(data.frame(t = res$statistic, pvalue = res$p.value))
})
hdac2 <- cbind(hdac2, tmp[, 11:12])
write.table(hdac2,  paste0(input_dir, "hdac2_DE.tsv"), quote = F, sep = "\t", row.names = F)

name_map <- read.csv("data/metabolomics/name_map.csv")
idx <- match(hdac2$Name, name_map$Query)
hdac2 <- cbind(hdac2, name_map[idx, c("Match", "HMDB", "PubChem", "KEGG")])

hdac2 <- hdac2 %>%
    mutate(Match = ifelse(is.na(Match), Name, Match)) %>%
    arrange(Match, pvalue) %>%
    filter(!duplicated(Match)) %>%
    mutate(regulate = ifelse(pvalue >= 0.05, "Normal", ifelse(log2FC > 0, "Up", "Down")), 
           log10P = -log10(pvalue))
write.table(hdac2,  paste0("res/3_epi/metabolomics/hdac2_DE_undup.tsv"), quote = F, sep = "\t", row.names = F)

mycol = brewer.pal(9, "Set1")[c(1, 2, 9)]
names(mycol) = c("Up", "Down", "Normal")
ggplot(hdac2, mapping = aes(x = log2FC, y = log10P, color = regulate)) + 
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = mycol) + 
    theme_bw() + 
    theme(text = tmp.text, axis.text.x = tmp.text, axis.text.y = tmp.text)
ggsave(file = paste0("plot/3_epi/metabolomics/", "HDAC2_H1299_Volcano.pdf"), width = 6, height = 5)


# Load PANC1_FOSL1 =======================================================
rm(list = ls())
tmp.text <- element_text(family="sans", size=12)
input_dir = "data/metabolomics/PANC1_FOSL1/"
header <- read.delim(paste0(input_dir, "header.txt"))
sam_cols <- c(paste0("QC", 1:3), paste0("SCR", 1:5), paste0("SH", 1:5))

pos <- read.delim(paste0(input_dir, "Area_0_202210211056.txt"), skip = 4, check.names = F, quote = "", comment.char = "") %>%
    mutate(ID = paste0(.data$`Average Mz`, "__", .data$`Average Rt(min)`*60), 
           ID2 = paste0("POS", .data$`Alignment ID`), 
           Name = sub(";.+", "", sub('"$', '', sub('^"', '', .data$`Metabolite name`)))) 
idx <- intersect(grep("Unknown", pos$Name, invert = T), grep("^w/o", pos$Name, invert = T))
pos <- pos[idx, c("ID2", "Name", sam_cols)]

neg <- read.delim(paste0(input_dir, "Area_0_20221021117.txt"), skip = 4, check.names = F, quote = "", comment.char = "") %>%
    mutate(ID = paste0(.data$`Average Mz`, "__", .data$`Average Rt(min)`*60), 
           ID2 = paste0("NEG", .data$`Alignment ID`),  
           Name = sub(";.+", "", sub('"$', '', sub('^"', '', .data$`Metabolite name`))))  
idx <- intersect(grep("Unknown", neg$Name, invert = T), grep("^w/o", neg$Name, invert = T))
neg <- neg[idx, c("ID2", "Name", sam_cols)]  

## DE
fosl1 <- rbind(pos, neg) 
fosl1$CV <- matrixStats::rowSds(as.matrix(fosl1[paste0("QC", 1:3)]))/
    rowMeans(as.matrix(fosl1[paste0("QC", 1:3)]))
fosl1 <- fosl1[!is.na(fosl1$CV) & fosl1$CV < 0.3, c("ID2", "Name", paste0("SCR", 1:5), paste0("SH", 1:5), "CV")]
fosl1_log <- log2(fosl1[, c(paste0("SCR", 1:5), paste0("SH", 1:5))] + 1)
fosl1$log2FC <- log2(rowMeans(fosl1[paste0("SH", 1:5)])/rowMeans(fosl1[paste0("SCR", 1:5)]))
tmp <- plyr::adply(fosl1_log, 1, function(x) {
    res = t.test(x[6:10], x[1:5])
    return(data.frame(t = res$statistic, pvalue = res$p.value))
})
fosl1 <- cbind(fosl1, tmp[, 11:12])
write.table(fosl1,  paste0(input_dir, "fosl1_DE.tsv"), quote = F, sep = "\t", row.names = F)

## map name 
name_map <- read.csv("data/metabolomics/name_map.csv")
idx <- match(fosl1$Name, name_map$Query)
fosl1 <- cbind(fosl1, name_map[idx, c("Match", "HMDB", "PubChem", "KEGG")])

fosl1 <- fosl1 %>%
    mutate(Match = ifelse(is.na(Match), Name, Match)) %>%
    arrange(Match, pvalue) %>%
    filter(!duplicated(Match)) %>%
    mutate(regulate = ifelse(pvalue >= 0.05, "Normal", ifelse(log2FC > 0, "Up", "Down")), 
           log10P = -log10(pvalue))
write.table(fosl1,  paste0("res/3_epi/metabolomics/fosl1_DE_undup.tsv"), quote = F, sep = "\t", row.names = F)

mycol = brewer.pal(9, "Set1")[c(1, 2, 9)]
names(mycol) = c("Up", "Down", "Normal")
ggplot(fosl1, mapping = aes(x = log2FC, y = log10P, color = regulate)) + 
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = mycol) + 
    theme_bw() + 
    theme(text = tmp.text, axis.text.x = tmp.text, axis.text.y = tmp.text)
ggsave(file = paste0("plot/3_epi/metabolomics/", "PANC1_FOSL1_Volcano.pdf"), width = 6, height = 5)



# compound <- c(hdac2$Name, fosl1$Name) %>% sort() %>% unique()
# write(compound, file = "data/metabolomics/all_compound.txt")
