# ==================================================================================================
# Author: Zhou Zhe
# Program: Generate MPs from Genes_nmf_w_basis, modified from Gavish et.al.
#          Multiple combinations of parameters were tested, and finally program size of 30 genes, 
#          intra-min of 21 genes (0.7), intra-max and inter-min of 6 genes (0.2) were chose. 
# Version: 1.0
# Date: Aug 8, 2023
# ==================================================================================================
    
# ----------------------------------------------------------------------------------------------------
# Load packages and functions  
# ----------------------------------------------------------------------------------------------------
rm(list = ls())
library(reshape2)
library(NMF)
library(ggplot2)
library(scales)
library(parallel)
library(clusterProfiler)
source("src/03_NMF/00_tools_custom_magma.R")
source("src/03_NMF/00_tools_robust_nmf_programs.R")   


## Input:
# Genes_nmf_w_basis is a list in which each entry contains NMF gene-scores of a single sample. In our study we ran NMF using ranks 4-9 on the top 7000 genes in each sample. Hence each entry in Genes_nmf_w_basis is a matrix with 7000 rows (genes) X 39 columns (NMF programs)  
# For the code below to run smoothly, please use the following nomenclature:
  # 1) End entry names in Genes_nmf_w_basis (i.e. each sample name) with '_rank4_9_nruns10.RDS' 
  # 2) End each matrix column name with an extension that represents the NMF rank and program index. for example '_rank4_9_nrun10.RDS.4.1' to represent the first NMF program obtained using rank=4, or '_rank4_9_nrun10.RDS.6.5' to represent the fifth NMF program obtained using rank=6    
# See Genes_nmf_w_basis_example.RDS for an example 
# We define MPs in 2 steps:
  # 1) The function robust_nmf_programs.R performs filtering, so that programs selected for defining MPs are:
  #    i) Robust                - recur in more that one rank within the sample 
  #    ii) Non-redundant        - once a NMF program is selected, other programs within the sample that are similar to it are removed
  #    iii) Not sample specific - has similarity to NMF programs in other samples 
  # ** Please see https://github.com/gabrielakinker/CCLE_heterogeneity for more details on how to define robust NMF programs 
  # 2) Selected NMFs are then clustered iteratively, as described in Figure S1. At the end of the process, each cluster generates a list of the 50 genes (i.e. the MP) that represent the NMF programs that contributed to the cluster. Notably, not all initially selected NMFs end up participating in a cluster  


# ----------------------------------------------------------------------------------------------------
# Select NMF programs
# ----------------------------------------------------------------------------------------------------


## Parameters 
all_celltypes <- c("Malignant", "Myeloid", "B_cells", "CD4_T", "CD8_T", "Mast_cells", "Endothelial", "Fibroblasts")
ct = "Malignant"
if (ct == "Malignant") {
    Genes_nmf_w_basis <- readRDS("res/02_NMF/NMF/malig/Genes_nmf_w_basis.RDS")
} else {
    Genes_nmf_w_basis <- readRDS(paste0("res/02_NMF/NMF/res_nonmalig/", ct, "_Genes_nmf_w_basis.RDS"))
}

program_size <- 30
intra_min_parameter <- 21
intra_max_parameter <- 6
inter_min_parameter <- 6

# get top 30 genes for each NMF program 
nmf_programs          <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:program_size]))
nmf_programs          <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 

# for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
nmf_filter_ccle       <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter=T, inter_min = inter_min_parameter)  
nmf_programs          <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs          <- do.call(cbind, nmf_programs)
nmf_programs_original <- nmf_programs

# calculate similarity between programs
cl = makeCluster(8, type = "FORK")
nmf_intersect         <- parApply(cl = cl, nmf_programs , 2, 
                                  function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y)))) 
stopCluster(cl)

# hierarchical clustering of the similarity matrix
nmf_intersect_hc     <- hclust(as.dist(program_size-nmf_intersect), method="average")
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]



# ----------------------------------------------------------------------------------------------------
# Cluster selected NMF programs to generate MPs
# ----------------------------------------------------------------------------------------------------


### Parameters for clustering
Min_intersect_initial <- 6    # the minimal intersection cutoff for defining the first NMF program in a cluster
Min_intersect_cluster <- 6    # the minimal intersection cutoff for adding a new NMF to the forming cluster 
Min_group_size        <- 5     # the minimal group size to consider for defining the first NMF program in a cluster 

Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list              <- list()   ### Every entry contains the NMFs of a chosen cluster
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()

nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1]>Min_group_size) {  
  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP                    <- nmf_programs[,names(Sorted_intersection[1])] # Genes in the forming MP are first chosen to be those in the first NMF. Genes_MP always has only 50 genes and evolves during the formation of the cluster
  nmf_programs                <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]  # remove selected NMF
  Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history                 <- Genes_MP  # has genes in all NMFs in the current cluster, for redefining Genes_MP after adding a new NMF 
  
  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.    
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {  
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[program_size])]   ### genes with overlap equal to the 50th gene
    
    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste( strsplit(i , "[.]")[[1]][1 : which(strsplit(i , "[.]")[[1]] == "RDS")]   , collapse = "."  )
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   
      
      Genes_MP_temp             <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[program_size])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:program_size] 
    }
    
    NMF_history     <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP        <- Genes_MP_temp[1:program_size]
    
    nmf_programs    <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]  # remove selected NMF
    
    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MMP_",k)]]           <- Genes_MP
  k <- k+1
  
  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # Remove current chosen cluster
  
  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   # Sort intersection of remaining NMFs not included in any of the previous clusters
  
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}


## Filter clusters ======
idx <- sapply(Cluster_list, length) >= 10
Cluster_list <- Cluster_list[idx]
names(Cluster_list) <- paste0("Cluster_", 1:length(Cluster_list))
MP_list <- MP_list[idx]
names(MP_list) <- paste0("MMP_", 1:length(MP_list))


## Re-Annotate MP genes by Zhou Zhe =======
library(scMetab)

TERM2GENE <- data.frame(TERM = gaudeMetabDf$Pathway, GENE = gaudeMetabDf$GeneSymbol)
MP_enrich <- lapply(MP_list, function(x) {
    en = enricher(x, TERM2GENE = TERM2GENE, minGSSize = 5, pvalueCutoff = 1, qvalueCutoff = 1)@result
    # en$geneRatio <- sapply(en$GeneRatio, function(text) eval(parse(text = text)))
    # en = en[order(en$geneRatio, decreasing = T), ]
})

msig <- rbind(read.gmt("/data/genome/MSigDB_GMT/msigdb_v2023.1.Hs_GMTs/c8.all.v2023.1.Hs.symbols.gmt"),
              read.gmt("/data/genome/MSigDB_GMT/msigdb_v2023.1.Hs_GMTs/c5.go.v2023.1.Hs.symbols.gmt"),
              read.gmt("/data/genome/MSigDB_GMT/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt"))
cl = makeCluster(10, type = "FORK")
MP_msig <- parLapply(cl = cl, MP_list, function(x) {
    en = enricher(x, TERM2GENE = msig, minGSSize = 5)@result
})
stopCluster(cl)

### get the source datasets of MMPs
Cluster_source <- lapply(Cluster_list, function(x) {
    x_split = strsplit(x, "_")
    res = sapply(x_split, function(y) {
        if (length(grep("valid", y)) > 0) paste(y[1:2], collapse = "_") else y[1]
    })
})
lapply(Cluster_source, table)


MP_intersection <- sapply(MP_list, function(x) sapply(MP_list, function(y) length(intersect(x, y))))


####  Sort Jaccard similarity plot according to new clusters:

inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))
  
}
inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered

nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new]) 

ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
ggsave(filename = paste0("plot/2.5_NMF/", ct, "_MP_cluster_p30_7_.2_final.pdf"), width = 10, height = 9)

MP_list <-  do.call(cbind, MP_list)
MP_list <- as.data.frame(MP_list)

save(Cluster_list, MP_list, Cluster_source, MP_enrich, MP_msig, MP_intersection, nmf_programs_original, nmf_intersect_original,
     file = paste0("res/2.5_NMF/MPs/", ct, "_MPs_p30_.7_.2_final.RData"))


