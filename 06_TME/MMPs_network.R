# ====================================================================
# Author: Zhou Zhe
# Function: MMP correlation of different cell types
# Version: 1.0
# Date: Aug 11, 2023
# ====================================================================

rm(list = ls())
library(plyr)
library(igraph)
library(RColorBrewer)
library(scMetab)
library(clusterProfiler)
library(parallel)

used.celltype = c("Malignant", "Myeloid", "B_cells", "CD4_T", "CD8_T", "Endothelial", "Fibroblasts")
sam_info_list <- readRDS("data/rds/NMF/cells_info_per_sam.RDS")
# celltype_color <- c(CD4_T = "#984EA3", CD8_T = "#984EA3", B_cells = "#FF7F00", celltype_color)
mycolor = c(scMetab::celltype_color[c(1, 3, 5)], brewer.pal(8, "Dark2")[3], 
            brewer.pal(12, "Paired")[10], scMetab::celltype_color[c(6, 7)])
names(mycolor) = used.celltype
MinScore = 0.5

mtx = lapply(used.celltype, function(ct) {
    
    MP_scores_per_sample <- readRDS(paste0("res/02_NMF/MPs/", ct, "_MP_scores_p30_.7_.2_final_MinScore_0.8.RDS"))
    adj_score <- plyr::ldply(MP_scores_per_sample, function(x) {
        colSums(x > MinScore)/nrow(x)
    })
    rownames(adj_score) <- adj_score$.id
    colnames(adj_score) <- c("samID", paste(ct, colnames(adj_score)[-1], sep = "_"))
    return(adj_score)
})
names(mtx) = used.celltype

#mtx <- data.frame(t(do.call(rbind, mtx)), check.names = F)
mtx <- Reduce(function(x, y) {
    merge(x, y, by = "samID", all = T)
}, mtx)
rownames(mtx) <- mtx$samID
mtx <- mtx[-1]

used_sam <- unique(do.call(rbind, lapply(sam_info_list, function(x) x$sam_info)))
used_sam <- used_sam[match(rownames(mtx), used_sam$samID), ]
all(rownames(mtx) == used_sam$samID)

mtx <- do.call(rbind, lapply(split(mtx, used_sam$Dataset), scale))


## correlation
ct_cor_pearson <- Hmisc::rcorr(mtx, type = "pearson")
ct_cor_spearman <- Hmisc::rcorr(mtx, type = "spearman")

cor_df = Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2"), all = T), 
                list(reshape2::melt(ct_cor_pearson$r), reshape2::melt(ct_cor_pearson$P), 
                     reshape2::melt(ct_cor_spearman$r), reshape2::melt(ct_cor_spearman$P)))
colnames(cor_df) = c("Vertice1", "Vertice2", "R", "P_pearson", "Rho", "P_spearman")

cor_df <- cor_df[!is.na(cor_df$P_pearson) & (cor_df$P_pearson < 0.05) & (cor_df$P_spearman < 0.05), ]
cor_df$logP = -log10(cor_df$P_pearson)
cor_df$ID <- apply(cor_df, 1, function(x) paste(sort(x[1:2]), collapse = "_"))
cor_df <- cor_df[!duplicated(cor_df$ID), ]


## plot for positive correlations 
pos = cor_df[(cor_df$R > 0) & (cor_df$logP > 5), c(1,2,7,3:6)]
pos_col_fun <- circlize::colorRamp2(c(4, 8, 12), brewer.pal(3, "YlOrRd"))
pos_ver = data.frame(Vertice = unique(c(pos$Vertice1, pos$Vertice2)))
pos_ver$celltype = sub("_MMP.+", "", pos_ver$Vertice)
pos_graph = graph_from_data_frame(pos, directed = F, vertices = pos_ver)
V(pos_graph)$color = mycolor[V(pos_graph)$celltype]
V(pos_graph)$label.color = "black"
V(pos_graph)$frame.width = 0
E(pos_graph)$width = 5
E(pos_graph)$color <- pos_col_fun(E(pos_graph)$logP)
pdf("plot/06_TME/TME_MMP_crosstalk_pos_sig5.pdf", width = 6, height = 6)
set.seed(2021)
V(pos_graph)$label = sub(".+_", "", names(V(pos_graph)))
plot(pos_graph, layout=layout.kamada.kawai, vertex.size=13, vertex.label.cex=0.5, 
     vertex.label.dist = 0)

V(pos_graph)$label = sub(" .+", "", V(pos_graph)$label)

plot.igraph(pos_graph, layout=layout.kamada.kawai, vertex.size=13, vertex.label.cex=0.5, 
     vertex.label.dist = 0)
dev.off()

library(ggplot2)
tmpdf = data.frame(x = 1:3, y = 1:3, value = c(4, 8, 12))
ggplot(tmpdf, aes(x = x, y = y, color = value)) + 
    geom_point() + 
    scale_color_gradientn(colours = brewer.pal(3, "YlOrRd"), breaks = c(4, 8, 12))
ggsave("plot/07_TME/TME_MMP_crosstalk_pos_sig5_legend.pdf")

## validate MMP clusters
MP_all <- readRDS("res/02_NMF/MPs/MP_list.RDS")
MP_cluster1 <- list(CD4_T_MMP2 = MP_all$CD4_T$`MMP2 Transport`, 
                    CD8_T_MMP2 = MP_all$CD8_T$`MMP2 Transport`, 
                    Myeloid_MMP2 = MP_all$Myeloid$MMP2, 
                    Endothelial_MMP3 = MP_all$Endothelial$MMP3, 
                    Fibroblasts_MMP4 = MP_all$Fibroblasts$MMP4)
jaccard_cluster1 <- sapply(MP_cluster1, function(x) {
    sapply(MP_cluster1, function(y) length(intersect(x, y))/length(union(x, y)))
})
MP_cluster1_pooled <- unique(unlist(MP_cluster1))

# Annotate MP genes ==============================
TERM2GENE <- data.frame(TERM = gaudeMetabDf$Pathway, GENE = gaudeMetabDf$GeneSymbol)
MP_enrich1 <- lapply(MP_cluster1, function(x) {
    en = enricher(x, TERM2GENE = TERM2GENE, minGSSize = 5, pvalueCutoff = 1, qvalueCutoff = 1)@result
    # en$geneRatio <- sapply(en$GeneRatio, function(text) eval(parse(text = text)))
    # en = en[order(en$geneRatio, decreasing = T), ]
})
MP_enrich1_pooled = enricher(unique(unlist(MP_cluster1)), TERM2GENE = TERM2GENE, minGSSize = 5, 
                             pvalueCutoff = 1, qvalueCutoff = 1)@result

msig <- rbind(#read.gmt("/mnt/c/data/genome/MSigDB_GMT/msigdb_v2023.1.Hs_GMTs/c8.all.v2023.1.Hs.symbols.gmt"),
              read.gmt("/mnt/c/data/genome/MSigDB_GMT/msigdb_v2023.1.Hs_GMTs/c5.go.bp.v2023.1.Hs.symbols.gmt"),
              read.gmt("/mnt/c/data/genome/MSigDB_GMT/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt"))
cl = makeCluster(10, type = "FORK")
MP_msig1 <- parLapply(cl = cl, MP_cluster1, function(x) {
    en = enricher(x, TERM2GENE = msig, minGSSize = 5)@result
})
stopCluster(cl) 

MP_msig1_pooled = enricher(unique(unlist(MP_cluster1)), 
                           TERM2GENE = msig, minGSSize = 5, qvalueCutoff = 0.01)@result
writexl::write_xlsx(MP_msig1_pooled, "res/07_TME/Pathway_enrichment_of_Cluster5.xlsx")


## plot for negative correlations
neg = cor_df[(cor_df$R < 0) & (cor_df$logP > 3), c(1,2,7,3:6)]
neg_col_fun <- circlize::colorRamp2(c(2:4), brewer.pal(3, "YlOrRd"))
neg_ver = data.frame(Vertice = unique(c(neg$Vertice1, neg$Vertice2)))
neg_ver$celltype = sub("_MMP.+", "", neg_ver$Vertice)
neg_graph = graph_from_data_frame(neg, directed = F, vertices = neg_ver)
V(neg_graph)$color = mycolor[V(neg_graph)$celltype]
V(neg_graph)$label.color = "black"
E(neg_graph)$width <- 5
V(neg_graph)$frame.width = 0
E(neg_graph)$color <- neg_col_fun(E(neg_graph)$logP)
pdf("plot/06_TME/TME_MMP_crosstalk_neg_sig3.pdf", width = 6, height = 6)
set.seed(2021)
V(neg_graph)$label = sub(".+_", "", names(V(neg_graph)))
plot(neg_graph, layout=layout.fruchterman.reingold, vertex.size=13, vertex.label.cex=0.5, 
     vertex.label.dist = 0)

V(neg_graph)$label = sub(" .+", "", V(neg_graph)$label)

plot.igraph(neg_graph, layout=layout.fruchterman.reingold, vertex.size=13, vertex.label.cex=0.5, 
            vertex.label.dist = 0)
dev.off()

