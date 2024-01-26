# ==================================================================================================
# Author: Zhou Zhe
# Program: Prepare expression matrix for NMF analysis, and process Genes_nmf_w_basis from NMF
# Version: 1.0
# Date: Jul 21, 2023
# ==================================================================================================

rm(list=ls())
library(Seurat)
library(scMetab)
library(NMF)
library(parallel)
# cl = makeCluster(8)


dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
non_malig <- c("Myeloid", "B cells", "CD4_T", "CD8_T", "Mast cells", "Endothelial", "Fibroblasts")

# I. Prepare metabolic gene expression matrix for NMF (malignant cells) ========
for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(Sys.time(), " ------ ", tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(paste0("data/rds/TpN/", tumor, ".TpN.rds"))
    metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu_obj))
    seu_obj <- seu_obj[metab.gene, ]
    
    if (!("samID" %in% colnames(seu_obj[[]]))) seu_obj$samID = seu_obj$patientID
    if (tumor == "PRAD") seu_obj$samID <- seu_obj$samID2
    
    epi <- subset(seu_obj, celltype == "Malignant" & TorN == "T")
    epi.list <- SplitObject(epi, split.by = "samID")
    nCells = sapply(epi.list, function(x) dim(x)[2])
    epi.list <- epi.list[nCells > 20]
    
    expr.list <- sapply(names(epi.list), function(sam) {
        print(paste0(Sys.time(), " ------ ", "Analyze sample ", sam))
        seu <- epi.list[[sam]]
        seu <- NormalizeData(seu)
        mean_expr <- sort(rowMeans(seu@assays$RNA@data), decreasing = T)
        kept_genes <- names(mean_expr)[1:500]
        seu <- ScaleData(seu, features = kept_genes)
        expr <- seu@assays$RNA@scale.data
        expr[expr < 0] = 0
        
        return(expr)
    })
    # names(expr.list) <- paste0(tumor, "_", names(epi.list), "_rank3_6_nruns10.RDS")
    saveRDS(expr.list, file = paste0("data/rds/NMF/malig_expr/", tumor, "_malig_expr_top500.rds"))
  
}


# II. Prepare metabolic gene expression matrix for NMF (non-malignant cells) ========
for (i in 1:nrow(dataset)) {
    tumor <- dataset$DataSets[i]
    print(paste0(Sys.time(), " ------ ", tumor, "========"))
    
    seu_obj <- readRDS(paste0("data/rds/TpN/", dataset$DataSets[i], ".TpN.rds"))
    metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu_obj))
    seu_obj <- seu_obj[metab.gene, seu_obj$TorN == "T"]
    
    non_malig <- intersect(c("Myeloid", "B cells", "CD4_T", "CD8_T", "Mast cells", "Endothelial", "Fibroblasts"), 
                           seu_obj$celltype2)
    non_malig.list <- lapply(non_malig, function(ct) {
        print(paste0(Sys.time(), " ------ ", "Analyze cell type: ", ct))
        ct_seu <- subset(seu_obj, celltype2 == ct)
        ct_list = SplitObject(ct_seu, split.by = "samID")
        nCells = sapply(ct_list, function(x) dim(x)[2])
        ct_list <- ct_list[nCells > 20]
        
        expr.list <- lapply(names(ct_list), function(sam) {
            print(paste0(" ------ ", Sys.time(), " ------ ", "Analyze sample ", sam))
            seu <- ct_list[[sam]]
            seu <- NormalizeData(seu, verbose = F)
            mean_expr <- sort(rowMeans(seu@assays$RNA@data), decreasing = T)
            kept_genes <- names(mean_expr)[1:min(500, sum(mean_expr > 0))]
            seu <- ScaleData(seu, features = kept_genes, verbose = F)
            expr <- seu@assays$RNA@scale.data
            expr[expr < 0] = 0
            
            return(expr)
        })
        names(expr.list) <- names(ct_list)
        return(expr.list)
    })
    names(non_malig.list) <- non_malig
    
    # names(expr.list) <- paste0(tumor, "_", names(epi.list), "_rank3_6_nruns10.RDS")
    saveRDS(non_malig.list, file = paste0("data/rds/NMF/non_malig/", tumor, "_nonmalig_expr_top500.rds"))
}


## save cell names and other info for each samples ------
malig_list <- lapply(dataset$DataSets, function(tumor) {
    readRDS(paste0("data/rds/NMF/malig_expr/", tumor, "_malig_expr_top500.rds"))
})
non_malig_list <- lapply(dataset$DataSets, function(tumor) {
    readRDS(paste0("data/rds/NMF/non_malig/", tumor, "_nonmalig_expr_top500.rds"))
})
names(non_malig_list) <- names(malig_list) <- dataset$DataSets

### save per cell type 
for(tumor in dataset$DataSets) non_malig_list[[tumor]]$Malignant = malig_list[[tumor]]
sam_info_list <- lapply(c("Malignant", non_malig), function(ct) {
    res <- lapply(dataset$DataSets, function(tumor) {
        x = non_malig_list[[tumor]]
        if (ct %in% names(x)) {
            x = x[[ct]]
            names(x) = paste0(tumor, "_", names(x))
            
            sam_info = data.frame(samID = names(x), 
                                  Tumor_Type = sub("_.+", "", tumor), 
                                  Dataset = tumor)
            sam_cells = lapply(x, function(y) colnames(y))
            return(list(sam_info = sam_info, sam_cells = sam_cells))
        } else return(list())
    })
    res1 = list(sam_info = do.call(rbind, lapply(res, function(x) x$sam_info)), 
                sam_cells = do.call(c, lapply(res, function(x) x$sam_cells)))
    return(res1)
})
names(sam_info_list) <- c("Malignant", non_malig)
saveRDS(sam_info_list, "data/rds/NMF/cells_info_per_sam.RDS")



# III. Prepare log2(cpm/10+1) matrix ========
all_cpm <- lapply(1:nrow(dataset), function(i) {
    tumor <- dataset$DataSets[i]
    print(paste0(Sys.time(), " ------ ", tumor, "========"))
    
    # I. Load Data ========
    seu_obj <- readRDS(paste0("data/rds/TpN/", dataset$DataSets[i], ".TpN.rds"))
    seu_obj <- seu_obj[, seu_obj$TorN == "T"]
    
    all_celltype <- intersect(c("Malignant", non_malig), seu_obj$celltype2)
    cpm.list <- lapply(all_celltype, function(ct) {
        print(paste0(Sys.time(), " ------ ", "Analyze cell type: ", ct))
        ct_seu <- subset(seu_obj, celltype2 == ct)
        ct_list = SplitObject(ct_seu, split.by = "samID")
        nCells = sapply(ct_list, function(x) dim(x)[2])
        ct_list <- ct_list[nCells > 20]
        
        expr.list <- lapply(names(ct_list), function(sam) {
            print(paste0(" ------ ", Sys.time(), " ------ ", "Analyze sample ", sam))
            seu <- ct_list[[sam]]
            seu <- NormalizeData(seu, scale.factor = 1e5, verbose = F)
            expr <- seu@assays$RNA@data / log(2)  ## transform to log2(CPM/10 + 1)
            return(expr)
        })
        names(expr.list) <- names(ct_list)
        return(expr.list)
    })
    names(cpm.list) <- all_celltype
    return(cpm.list)
})
names(all_cpm) <- dataset$DataSets

### save per cell type 
cpm_per_celltype <- lapply(c("Malignant", non_malig), function(ct) {
    res <- lapply(dataset$DataSets, function(tumor) {
        x = all_cpm[[tumor]]
        if (ct %in% names(x)) {
            x = x[[ct]]
            names(x) = paste0(tumor, "_", names(x))
            return(x)
        } else return(list())
    })
    res <- do.call(c, res)
    fastSave::saveRDS.pigz(res, paste0("data/rds/NMF/log_transformed/", sub(" ", "_", ct), "_log2_all_samples.RDS"), n.cores = 10)
    return(res)
})
names(cpm_per_celltype) <- c("Malignant", non_malig)


# IV. run NMF on server ========
## see 00_tools_run_NMF.R

# V. Prepare Genes_nmf_w_basis for malignant cells ========
f_sam <- dir("res/2.5_NMF/NMF/malig/", pattern = "_rank3_6_nruns10.RDS")
Genes_nmf_w_basis_rank3_6 <- lapply(f_sam, function(sam) {
    nmf_res <- readRDS(file.path("res/2.5_NMF/malig/", sam))
    w_basis = sapply(3:6, function(n_rank) {
        w = NMF::basis(nmf_res$fit[[as.character(n_rank)]])
        colnames(w) = paste(sam, n_rank, 1:n_rank, sep = ".")
        return(w)
    })
    w_basis = do.call(cbind, w_basis)
})
names(Genes_nmf_w_basis_rank3_6) <- f_sam
saveRDS(Genes_nmf_w_basis_rank3_6, "res/2.5_NMF/NMF/malig/Genes_nmf_w_basis_rank3_6.RDS")

f_sam <- dir("res/2.5_NMF/NMF/malig/", pattern = "_rank7_9_nruns10.RDS")
Genes_nmf_w_basis_rank7_9 <- lapply(f_sam, function(sam) {
    nmf_res <- readRDS(file.path("res/2.5_NMF/malig/", sam))
    w_basis = sapply(7:9, function(n_rank) {
        w = NMF::basis(nmf_res$fit[[as.character(n_rank)]])
        colnames(w) = paste(sam, n_rank, 1:n_rank, sep = ".")
        return(w)
    })
    w_basis = do.call(cbind, w_basis)
})
names(Genes_nmf_w_basis_rank7_9) <- f_sam
saveRDS(Genes_nmf_w_basis_rank7_9, "res/2.5_NMF/NMF/malig/Genes_nmf_w_basis_rank7_9.RDS")

all_sams = sub("_rank7_9_nruns10.RDS", "", f_sam)
names(Genes_nmf_w_basis_rank3_6) <- sub("_rank3_6.+", "", names(Genes_nmf_w_basis_rank3_6))
names(Genes_nmf_w_basis_rank7_9) <- sub("_rank7_9.+", "", names(Genes_nmf_w_basis_rank7_9))
all(names(Genes_nmf_w_basis_rank3_6) == names(Genes_nmf_w_basis_rank7_9))

Genes_nmf_w_basis <- lapply(all_sams, function(sam) {
    a = Genes_nmf_w_basis_rank3_6[[sam]]
    b = Genes_nmf_w_basis_rank7_9[[sam]]
    if (all(rownames(a) == rownames(b))) res = cbind(a, b) else stop("Genes not the same!")
    res <- res[, -c(1:3)]
    colnames(res) <- sub("rank3_6", "rank4_9", x = colnames(res))
    colnames(res) <- sub("rank7_9", "rank4_9", x = colnames(res))
    return(res)
})
names(Genes_nmf_w_basis) <- paste0(all_sams, "_rank4_9_nruns10.RDS")
saveRDS(Genes_nmf_w_basis, "res/2.5_NMF/NMF/malig/Genes_nmf_w_basis.RDS")
# stopCluster(cl)


# VI. Prepare Genes_nmf_w_basis for non-malignant cells ========
Genes_nmf_w_basis_list <- lapply(non_malig, function(ct) {
    ct <- sub(" ", "_", ct)
    f_sam <- dir("res/2.5_NMF/NMF/res_nonmalig/", pattern = paste0(ct, "_rank4_9_nruns10.RDS"))
    Genes_nmf_w_basis <- lapply(f_sam, function(sam) {
        nmf_res <- readRDS(file.path("res/2.5_NMF/NMF/res_nonmalig/", sam))
        w_basis = sapply(4:9, function(n_rank) {
            w = NMF::basis(nmf_res$fit[[as.character(n_rank)]])
            colnames(w) = paste(sam, n_rank, 1:n_rank, sep = ".")
            return(w)
        })
        w_basis = do.call(cbind, w_basis)
    })
    names(Genes_nmf_w_basis) <- f_sam
    saveRDS(Genes_nmf_w_basis, paste0("res/2.5_NMF/NMF/res_nonmalig/", ct, "_Genes_nmf_w_basis.RDS"))
    return(Genes_nmf_w_basis)
})
names(Genes_nmf_w_basis_list) <- non_malig 
