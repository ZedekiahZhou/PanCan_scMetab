# ==================================================================================================
# Author: Zhou Zhe
# Program: Run NMF (on server)
# Version: 1.0
# Date: Jul 21, 2023
# ==================================================================================================

library(NMF)

used = commandArgs(trailingOnly = T)
dt = list(datasets1 = c('BRCA', 'BRCA_valid2', 'LUAD_valid2', 'PAAD'), 
          datasets2 = c('BRCA_valid1', 'LUAD_valid1'),
          datasets3 = c('COAD', 'LUAD'),
          datasets4 = c('PRAD', 'STAD'))
datasets = dt[[used]]

fpath = "data/non_malig/"

# sink(paste0("./", used, ".log"))
for (tumor in datasets) {
    print(paste0(Sys.time(), " ------ ", tumor, " ========"))
    
    expr_list = readRDS(file = paste0(fpath, tumor, "_nonmalig_expr_top500.rds"))
    mnf_list = lapply(names(expr_list), function(ct) {
        print(paste0(" ====== ", Sys.time(), " ------ ", ct, " --------"))
        ct_list = expr_list[[ct]]
        res = sapply(names(ct_list), function(sam) {
            print(paste0("            ", Sys.time(), " ------ ", sam, " --------"))
            expr = ct_list[[sam]]
            nmf_res <- nmf(expr, rank = 4:9, nrun = 10, method = "snmf/r")
            saveRDS(nmf_res, file = paste0("res_nonmalig/", tumor, "_", sam, "_", sub(" ", "_", ct), 
                                           "_rank4_9_nruns10.RDS"))
        })
    })
        
        
}
# sink()
