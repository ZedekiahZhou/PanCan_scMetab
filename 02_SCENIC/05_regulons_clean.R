# ====================================================================
# Author: Zhou Zhe
# Function: Regulons data
# Version: 1.0
# Date: Jul 11, 2022
# ====================================================================

rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(patchwork)
library(GSEABase)
library(clusterProfiler)
library(ggvenn)
detach("package:scMetab", unload = TRUE)
library(scMetab)
library(parallel)
library(networkD3)

pseudo <- 1e-50
mycol <- brewer.pal(5, "RdBu")

dir.used <- "5_SCENIC/"
file.prefix <- "Regulon_"

dataset <- readxl::read_excel("../Datasets.xlsx")

# function to calculate jaccard index of two gene sets list
getJacc <- function(gs1, gs2) {
    gs1 = geneIds(gs1)
    gs2 = geneIds(gs2)
    common_tf = intersect(names(gs1), names(gs2))
    
    res = sapply(common_tf, function(x) jaccard_dist(gs1[[x]], gs2[[x]]))
    names(res) = common_tf
    return(res)
}

# metabolic genes
metab.genes <- unique(gaudeMetabDf$GeneSymbol)
metab.list <- split(gaudeMetabDf$GeneSymbol, f = gaudeMetabDf$Pathway)
metab.len <- sapply(metab.list, length)
metab.list <- metab.list[metab.len >=5]

metab.gmt <- lapply(names(metab.list), function(x) {
    GeneSet(metab.list[[x]], setName = x)
})
metab.gmt <- GeneSetCollection(metab.gmt)
toGmt(metab.gmt, "data/GeneSets/gaudeMetab.gmt")

# regulon list ======
regulon.list <- mclapply(1:nrow(dataset), function(i) {
    tumor = dataset$DataSets[i]
    regulons <- getGmt(paste0("res/02_SCENIC/regulon/", tumor, "_reg.gmt"))
    # regulons <- geneIds(regulons)
    
    reg_names <- sapply(regulons, setName)
    reg_names <- sub("(+)", "", reg_names, fixed = T)
    
    # filter out metabolic genes
    regulons_metab <- GeneSetCollection(lapply(1:length(reg_names), function(i) {
        x = regulons[[i]] & metab.genes
        setName(x) = reg_names[i]
        x
    }))
    
    # filter by length
    reg_len <- sapply(regulons_metab, function(x) length(geneIds(x)))
    regulons_metab <- regulons_metab[reg_len >= 5]
    # toGmt(regulons_metab, paste0("res/5_SCENIC/regulon/", tumor, "_reg_metab.gmt"))

    return(regulons_metab)
}, mc.cores = 1)
names(regulon.list) <- dataset$DataSets


# check pathway enrichment for each TF
TERM2GENE <- data.frame(TERM = gaudeMetabDf$Pathway, GENE = gaudeMetabDf$GeneSymbol)
TERM2GENE <- TERM2GENE[TERM2GENE$TERM %in% names(metab.gmt), ]
path_enrich.list <- mclapply(regulon.list, function(regulons) {
    res = lapply(regulons, function(x) {
        en = enricher(gene = geneIds(x), TERM2GENE = TERM2GENE, minGSSize = 5, pvalueCutoff = 1, qvalueCutoff = 1)@result
        en$geneRatio <- sapply(en$GeneRatio, function(text) eval(parse(text = text)))
        en = en[order(en$geneRatio, decreasing = T), ]
        # en = en[en$geneRatio > 0.1, ]
    })
    names(res) = names(regulons)
    return(res)
}, mc.cores = 1)

tmp <- lapply(path_enrich.list, function(x) x$FOSL1)

fosl1 <- path_enrich.list$LUAD$FOSL1
hdac2 <- path_enrich.list$LUAD$HDAC2

# fosl1 <- subset(fosl1, Count > 1)
fosl1 <- fosl1[1:15, ]
fosl1_plot <- data.frame(from = "FOSL1_reg", to = fosl1$ID, value = fosl1$Count)
fosl1_plot <- rbind(fosl1_plot, data.frame(from = "FOSL1_reg", to = "Other Pathways", 
                                           value = 53 - length(unique(unlist(strsplit(fosl1$geneID, "/", fixed = TRUE))))))
fosl1_nodes <- data.frame(ID = c("FOSL1_reg", fosl1_plot$to))
fosl1_plot$IDfrom <- match(fosl1_plot$from, fosl1_nodes$ID) - 1
fosl1_plot$IDto <- match(fosl1_plot$to, fosl1_nodes$ID) - 1
sankeyNetwork(Links = fosl1_plot, Nodes = fosl1_nodes, 
              Source = "IDfrom", Target = "IDto", Value = "value", NodeID = "ID", 
              fontFamily = "sans", fontSize = 11, height = 96*3, width = 96*5)

# hdac2 <- subset(hdac2, Count > 1)
hdac2 <- hdac2[1:15, ]
hdac2_plot <- data.frame(from = "HDAC2_reg", to = hdac2$ID, value = hdac2$Count)
hdac2_plot <- rbind(hdac2_plot, data.frame(from = "HDAC2_reg", to = "Other Pathways", 
                                           value = 116 - length(unique(unlist(strsplit(hdac2$geneID, "/", fixed = TRUE))))))
hdac2_nodes <- data.frame(ID = c("HDAC2_reg", hdac2_plot$to))
hdac2_plot$IDfrom <- match(hdac2_plot$from, hdac2_nodes$ID) - 1
hdac2_plot$IDto <- match(hdac2_plot$to, hdac2_nodes$ID) - 1
sankeyNetwork(Links = hdac2_plot, Nodes = hdac2_nodes, 
              Source = "IDfrom", Target = "IDto", Value = "value", NodeID = "ID", 
              fontFamily = "sans", fontSize = 11, height = 96*3, width = 96*5)


# check TF enrichment for each pathway
path_TF.list <- mclapply(regulon.list, function(regulons) {
    regulons_gene <- lapply(regulons, geneIds)
    names(regulons_gene) <- names(regulons)
    TERM2GENE <- reshape2::melt(regulons_gene)[c(2:1)]
    colnames(TERM2GENE) = c("TERM", "GENE")
    
    res = lapply(metab.list, function(x) {
        en = enricher(gene = x, TERM2GENE = TERM2GENE, minGSSize = 5)@result
        en$geneRatio <- sapply(en$GeneRatio, function(text) eval(parse(text = text)))
        en = en[order(en$geneRatio, decreasing = T), ]
        en = en[en$geneRatio > 0.1, ]
    })
}, mc.cores = 10)


# save enrich result
fastSave::save.pigz(path_enrich.list, path_TF.list, file = "res/5_SCENIC/pathway_TF_enrichment.rda")

# tmp <- lapply(path_TF.list, function(x) x$`Glycolysis and Gluconeogenesis`)
# sort(table(unlist(lapply(tmp, rownames))), decreasing = T)
# tmp2 <- lapply(path_TF.list, function(x) x$`Oxidative Phosphorylation`)
# sort(table(unlist(lapply(tmp2, rownames))), decreasing = T)
# 
# tmp1 <- lapply(regulon.list, function(x) geneIds(x[["MXI1"]]))




# TF.list <- lapply(regulon.list, names)
# shared.TF <- Reduce(intersect, TF.list)
# shared.brca <- Reduce(intersect, TF.list[1:3])
# shared.luad <- Reduce(intersect, TF.list[4:6])
# 
# tmpfile <- paste0("plot/", dir.used, file.prefix, "TF_venn.pdf")
# pdf(tmpfile, width = 3, height = 3)
# ggvenn(TF.list[c(1, 4, 7, 8, 9)], fill_color = brewer.pal(4, "Set2"), show_percentage = FALSE, 
#        set_name_size = 2.82, text_size = 2.82)
# dev.off()
