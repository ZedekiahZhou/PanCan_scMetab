# ==================================================================================================
# Author: Zhou Zhe
# Program: Program for depicting the MP distribution of a given cell types across and within samples 
#          (Modified from Gavish et. al.)
# Version: 1.0
# Date: Aug 11, 2023
# ==================================================================================================



# -------------------------------------------------------------------------------------------
# Program for depicting the MP distribution of a given cell types across and within samples 
# ------------------------------------------------------------------------------------------- 

# My_study  = a list in which each element represents a sample in the study; the elements are expression matrices of a given cell type (i.e. cell_type below) in CPM or TPM units. Row names in each matrix should be gene symbols, and column names the cell IDs. 
# cell_type = one of the following cell types : "Cancer" , "Endothelial" , "Epithelial" , "Fibroblasts" , "Macrophages" , "CD4_T_cells", "CD8_T_cells" , "B_cells". Here "Cancer" represents malignant cells
# MinGenes  = the min number of MP genes that should exist in the study (default is 25)
# MinScore  = the minimal score for assigning a cell to a MP (default is 1). A cell is assigned to the MP with the maximal score, given that it exceeded MinScore. If the maximal score is below MinScore, the cell is unassigned. 
# MinCells  = the minimal % of cells that should be assigned to a MP in order to account for the MP (default is 0.05)
# MP_list   = a list of the MPs for each cell_type

### Output: 
# (a) Pie chart of the MP distribution across the study
# (b) Bar plot per tumor with the % of cells per MP
# (c) Expression heatmap per tumor 



# ------------------------------------------------------------------------------------------- 

rm(list = ls())
library(ggplot2)
library(ggforce)
library(scalop)  # see https://rdrr.io/github/jlaffy/scalop/man/sigScores.html
library(gridExtra)
library(ggpubr)
library(scales)
library(Matrix)
library(parallel)



### Define the following:
celltype = "Malignant"
My_study_sparse  = readRDS(paste0("data/rds/NMF/log_transformed/", sub(" ", "_", celltype), "_log2_all_samples.RDS"))
# My_study = lapply(My_study_sparse, as.matrix)
para = commandArgs(trailingOnly = T)[1]  # such as "p30_.7_.2"
MinScore = as.numeric(commandArgs(trailingOnly = T)[2]) # c(0.5, 0.8, 1)

para = "p30_.7_.2_final"
MinScore = 0.8
load(paste0("res/02_NMF/MPs/", celltype, "_MPs_", para, ".RData"))
MP_all <- readRDS("res/02_NMF/MPs/MP_list.RDS")
heatCols  = readRDS("src/02_NMF/heatCols.RDS.gz")  # can be found in the Github repository


Score_cells <- function(L,
                        cell_type,
                        MinGenes = 25,
                        MP_list
) 
{
    
    MP_list <- MP_list[[cell_type]]
    
    withr::local_seed(2021)
    cl = makeForkCluster(nnodes = 4)
    system.time(
    MP_scores_per_sample    <- parLapply(cl, L, function(x)  {
        x = as.matrix(x)
        res = scalop::sigScores(x, sigs = MP_list, conserved.genes = MinGenes/50)
        gc()
        return(res)
    })
    )
    stopCluster(cl)
    return(MP_scores_per_sample)
}

Score_cells_win <- function(L,
                        cell_type,
                        MinGenes = 25,
                        MP_list
) 
{
    
    MP_list <- MP_list[[cell_type]]
    
    withr::local_seed(2021)
    cl = makeCluster(4)
    system.time(
        MP_scores_per_sample    <- parLapply(cl, L, function(x)  {
            x = as.matrix(x)
            res = scalop::sigScores(x, sigs = MP_list, conserved.genes = MinGenes/50)
            gc()
            return(res)
        })
    )
    stopCluster(cl)
    return(MP_scores_per_sample)
}


filter_score <- function(MP_scores_per_sample, 
                         MinScore = 1,
                         MinCells = 0.05
) {
    remove_cells <- function(x,MinScore){         # remove cells whose max score was below MinScore
        max_score <- apply(x, 1, function(y) max(y))
        cells_rm  <- which(max_score < MinScore)
        if (length(cells_rm) > 0 ){
            x <- x[-cells_rm , ]
        }
        return(x)
    }
    MP_scores_per_sample <- lapply(MP_scores_per_sample, function(x) remove_cells(x,MinScore))
    
    Assign_MP_per_cell   <- lapply(MP_scores_per_sample, function(x) apply(x, 1, function(y) colnames(x)[which(y==max(y))] ) ) 
    
    filter_cells <- function(x,MinCells){         # remove MP that were assassin to less than MinCells in each sample
        MP_frequency <- as.numeric(ave(x, x, FUN = length))
        MP_rm        <- which(MP_frequency/length(x) < MinCells)  # MPs to be removed
        if (length(MP_rm)>0){
            x <- x[-MP_rm]
        }
        return(x)
    }
    Assign_MP_per_cell_filtered <- lapply(Assign_MP_per_cell, function(x) filter_cells(x,MinCells)) 
    
    return(Assign_MP_per_cell_filtered)
}        
    

system.time(MP_scores_per_sample <- Score_cells_win(L = My_study_sparse , cell_type = celltype , MP_list=MP_all))
saveRDS(MP_scores_per_sample, file = paste0("res/02_NMF/MPs/", celltype, "_MP_scores_", para, "_MinScore_", MinScore, ".RDS"))

Assign_MP_per_cell_filtered <- filter_score(MP_scores_per_sample, MinScore = MinScore)
Assign_MP_per_cell_filtered <- Assign_MP_per_cell_filtered[sapply(Assign_MP_per_cell_filtered, length) > 0]
# saveRDS(Assign_MP_per_cell_filtered, paste0("res/2.5_NMF/MPs/Malignant_MP_Assigned", para, ".RDS"))



#### Part A: pie chart of MPs across the whole study
MinCells <- 0.05 
MP_list <- MP_all[[celltype]]
MP_color = c(colorRampPalette(scMetab::celltype_color)(length(MP_list)), "grey50")
names(MP_color) <- c(names(MP_list), "Other MPs (<5%)")

tmp <- table(unlist(Assign_MP_per_cell_filtered))
df       <- data.frame(MPs = names(tmp)  ,  frequency = as.vector(tmp)/sum(tmp))
MPs_rm   <- which(df$frequency/sum(df$frequency) < MinCells)  ## remove also MPs that were assassin to less than MinCells from the total  
if (length(MPs_rm > 0)){
    df$MPs[MPs_rm] <- paste0("Other MPs (<", MinCells*100, "%)")
}
df <- df[order(df$frequency, decreasing = T), ]
df$MPs <- factor(df$MPs, levels = unique(df$MPs))



ggplot(df) + 
    geom_arc_bar(stat = "pie", aes(x0=0,y0=0,r0=4,r=5, amount = frequency, fill = MPs), color = "white") + 
    scale_fill_manual(values = MP_color[names(MP_color) %in% df$MPs]) + 
    theme_void()
ggsave(paste0("plot/2.5_NMF/", celltype, "_MP_Assigned_pie_", para, "_MinScore_", MinScore, ".pdf"), width = 6, height = 4)

#### Part A1: pie chart of MPs each study
sam_info_list <- readRDS("data/rds/NMF/cells_info_per_sam.RDS")[[celltype]]
attach(sam_info_list)
My_study_sparse <- My_study_sparse[names(Assign_MP_per_cell_filtered)]
sam_info <- sam_info[match(names(Assign_MP_per_cell_filtered), sam_info$samID), ]


cancer_list <- tapply(Assign_MP_per_cell_filtered, sam_info$Tumor_Type, function(x) {
    y = table(unlist(x))
    tmp_df <- data.frame(MPs = names(y), frequency = as.vector(y)/sum(y))
    tmp_df <- tmp_df[order(tmp_df$frequency, decreasing = T), ]
    MPs_rm   <- which(tmp_df$frequency < MinCells)  ## remove also MPs that were assassin to less than MinCells from the total  
    if (length(MPs_rm > 0)){
        tmp_df$MPs[MPs_rm] <- paste0("Other MPs (<", MinCells*100, "%)")
    }
    return(tmp_df)
})

pdf(paste0("plot/2.5_NMF/", celltype, "_MP_Assigned_pie_tumor_type.pdf"), width = 6, height = 4.3)
cancer_plot <- lapply(names(cancer_list), function(cancer) {
    x = cancer_list[[cancer]]
    tmp_plot = ggplot(x) + 
        geom_arc_bar(stat = "pie", aes(x0=0,y0=0,r0=4,r=5, amount = frequency, fill = MPs), color = "white") + 
        scale_fill_manual(values = MP_color[names(MP_color) %in% x$MPs]) + 
        theme_void() + 
        ggtitle(cancer)
    print(tmp_plot)
    return(tmp_plot)
})
dev.off()

dataset_list <- tapply(Assign_MP_per_cell_filtered, sam_info$Dataset, function(x) {
    y = table(unlist(x))
    tmp_df <- data.frame(MPs = names(y), frequency = as.vector(y)/sum(y))
    tmp_df <- tmp_df[order(tmp_df$frequency, decreasing = T), ]
    MPs_rm   <- which(tmp_df$frequency < MinCells)  ## remove also MPs that were assassin to less than MinCells from the total  
    if (length(MPs_rm > 0)){
        tmp_df$MPs[MPs_rm] <- paste0("Other MPs (<", MinCells*100, "%)")
    }
    return(tmp_df)
})

pdf(paste0("plot/2.5_NMF/", celltype, "_MP_Assigned_pie_dataset.pdf"), width = 6, height = 4.3)
dataset_plot <- lapply(names(dataset_list), function(dataset) {
    x = dataset_list[[dataset]]
    tmp_plot = ggplot(x) + 
        geom_arc_bar(stat = "pie", aes(x0=0,y0=0,r0=4,r=5, amount = frequency, fill = MPs), color = "white") + 
        scale_fill_manual(values = MP_color[names(MP_color) %in% x$MPs]) + 
        theme_void() + 
        ggtitle(dataset)
    print(tmp_plot)
    return(tmp_plot)
})
dev.off()
    
#### Part B: Bar plot with MP distribution per sample

sample_num <- length(Assign_MP_per_cell_filtered)
nCol       <- floor(sqrt(sample_num))

### For arranging the full MP names as legend in the plot
MPs <- unique(unlist(Assign_MP_per_cell_filtered))    

Assign_MP_per_cell_filtered_abbrev <- lapply(Assign_MP_per_cell_filtered, function(x) apply(as.data.frame(x), 1, function(x)  strsplit(x, "[ ]")[[1]][1]))

df_list    <- lapply(Assign_MP_per_cell_filtered_abbrev, function(x) data.frame(MPs = names(table(x))  ,  frequency = as.vector(table(x))/sum(as.vector(table(x))) ))
df_list2   <- lapply(Assign_MP_per_cell_filtered, function(x) data.frame(MPs = names(table(x))  ,  frequency = as.vector(table(x))/sum(as.vector(table(x))) ))

P <- lapply(seq_along(df_list), function(I) 
    ggplot(data=df_list[[I]], aes(x=reorder(MPs,-frequency), y=frequency)) +
        geom_bar(stat="identity")+
        theme(axis.title.x=element_blank())+
        ggtitle(names(df_list)[I])
)
names(P) = names(df_list2)

# P

#### Part C: an expression heatmap for each sample

df_list_ordered  <- lapply(df_list, function(x) x[order(x$frequency,decreasing = T),])  ## order MPs in each sample as in the bar plots, in decreasing order
df_list_ordered2 <- lapply(df_list2, function(x) x[order(x$frequency,decreasing = T),])  ## order MPs in each sample as in the bar plots, in decreasing order
 
### extract and sort cells
L1        <- lapply(seq_along(My_study_sparse), function(I) My_study_sparse[[I]][ , names(Assign_MP_per_cell_filtered[[I]]), drop = F]) ## extract relevant cells that were assigned to an MP (at least 5%)
names(L1) <- names(My_study_sparse)
L2 <- lapply(seq_along(L1), function(I) L1[[I]][ , order(match(Assign_MP_per_cell_filtered_abbrev[[I]],df_list_ordered[[I]]$MPs)), drop = F] )  # sort cells according to MPs in df_list_ordered
names(L2) <- names(L1)

### select MP genes matrix and centering
MP_genes_per_sample <- lapply(df_list_ordered2, function(x) unlist(MP_list[x$MPs]))  ### MP genes sorted by the MP order in df_list_ordered2
L_plot        <- lapply(seq_along(L2), function(I) {
    idx = match(MP_genes_per_sample[[I]] , rownames(L2[[I]]))
    res = as.matrix(L2[[I]][idx[!is.na(idx)] , ])
    res = res - rowMeans(res)
}) 
names(L_plot) <- names(L2)

### plot function
color.scheme <- colorRampPalette(c(heatCols))(n=333)
Plot_heatmap <- function(M){ 
    
    M_new2        <- M
    M_new2        <- apply(M_new2, 2, rev)
    M_meltII      <-  reshape2::melt(t(M_new2)) 
    M_meltII$Var2 <- factor(M_meltII$Var2)
    
    G <- ggplot(data = M_meltII, aes(x=Var1, y=Var2, fill=value, color=value)) + 
        geom_raster() + 
        scale_color_gradient2(limits=c(-4,4), low=color.scheme[1:111],  mid =color.scheme[112:222], high = color.scheme[223:333], midpoint = 0, oob=squish, name=NULL) +                                
        scale_fill_gradient2(limits=c(-4,4), low=color.scheme[1:111],  mid =color.scheme[112:222], high = color.scheme[223:333], midpoint = 0, oob=squish, name=NULL)  +
        theme(  panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 8), legend.title = element_text(size=8), legend.text = element_text(size = 8), legend.text.align = 0.5, legend.justification = "bottom" ) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank() )
    
    return(G)
}

P1 <- lapply(L_plot, function(x) Plot_heatmap(x))

### combine bar plot and heatmap
P_plot <- c()
for (i in 1:length(P)){
    P_plot <- c(P_plot , P[i] , P1[i])
}

pdf(paste0("plot/2.5_NMF/", celltype, "_Sample_MPs_heatmap_raster_", para, "_MinScore_", MinScore, ".pdf"), width = 16, height = 9)
P_plot <- do.call("ggarrange", c(P_plot, ncol=2, nrow = 3))
lapply(P_plot, function(x) 
    
    print(annotate_figure(x, left   = text_grob(paste(MPs , collapse = "\n"), 
                                          color = "red", face = "bold", size = 10)))
    
)
dev.off()


