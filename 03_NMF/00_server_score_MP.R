library(parallel)


Score_cells <- function(L,
                        cell_type,
                        MinGenes = 25,
                        MP_list, 
                        nnodes = 10
) 
{
    
    MP_list <- MP_list[[cell_type]]
    
    withr::local_seed(2021)
    cl = makeForkCluster(nnodes = nnodes)
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

celltype = "Malignant"
My_study_sparse  = readRDS(paste0("data/rds/NMF/log_transformed/", sub(" ", "_", celltype), "_log2_all_samples.RDS"))

MP_all <- readRDS("res/02_NMF/MPs/MP_list.RDS")
MP_Gavish <- readRDS("res/02_NMF/MPs/Gavish_MP_list.RDS")


system.time(MP_scores_per_sample <- Score_cells(L = My_study_sparse , cell_type = celltype , MP_list=MP_all))
saveRDS(MP_scores_per_sample, file = paste0("res/02_NMF/MPs/", celltype, "_MP_scores_", para, "_MinScore_", MinScore, "_1107.RDS"))
