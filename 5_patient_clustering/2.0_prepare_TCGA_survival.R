# ==================================================================================================
# Author: Zhou Zhe
# Program: Clean and prepare gene expression and clinical data for TCGA survival analysis
# Version: 1.0
# Date: Aug 9, 2022
# ==================================================================================================

rm(list=ls())
library(SummarizedExperiment)
library(data.table)
library(scMetab)
library(GSEABase)
library(parallel)


# I. Prepare data ========
dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
dataset <- unique(dataset[, c("Tissue", "tumorType")])
dir.used <- "6_patient_clustering/SurvPre/"

# TCGA sample phenotype
pheno <- read.delim(gzfile(paste0("data/xena_tpm/TcgaTargetGTEX_phenotype.txt.gz")))
pheno$patient <- sub("-[0-9]+$", "", pheno$sample)
pheno <- pheno[pheno$X_sample_type == "Primary Tumor", ]  # No duplicate patients

# clinical data 
clinical_raw <- as.data.frame(readxl::read_excel("data/TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR"))
sam_clin <- merge(pheno, clinical_raw, by.x = "patient", by.y = "bcr_patient_barcode", all = FALSE)
rownames(sam_clin) <- sam_clin$patient

# ID to Symbol
id2symbol <- read.delim("data/xena_tpm/probeMap_gencode.v23.annotation.gene.probemap")
id2symbol$updated <- update_symbols(id2symbol$gene)$updated
id2symbol <- id2symbol[!duplicated(id2symbol$updated), ]

# do not use PRAD
all_clin <- lapply(1:5, function(i) {
    tumor = dataset$tumorType[i]
    tpm_data <- data.frame(fread(paste0("data/xena_tpm/", dataset$Tissue[i], "_matrix.tsv.gz")), 
                           row.names = 1, check.names = F)
    
    # only keep samples with both expression data and clinical data ---
    tumor_sam <- intersect(colnames(tpm_data), sam_clin$sample)
    tpm_data <- tpm_data[, tumor_sam]
    colnames(tpm_data) <- sam_clin$patient[match(tumor_sam, sam_clin$sample)]
    tumor_clin <- sam_clin[colnames(tpm_data), ]
    
    # filter gene --
    tpm_data <- tpm_data[rownames(tpm_data) %in% id2symbol$id, ]
    rownames(tpm_data) <- id2symbol$updated[match(rownames(tpm_data), id2symbol$id)]
    
    
    # filter by other clinical info ------
    tumor_clin <- tumor_clin[!is.na(tumor_clin$OS.time), ]
    tumor_clin <- tumor_clin[!is.na(tumor_clin$age_at_initial_pathologic_diagnosis), ]
    tumor_clin$age = tumor_clin$age_at_initial_pathologic_diagnosis
    tumor_clin$cancerType = tumor_clin$type
    tumor_clin$OS.time.month <- round(tumor_clin$OS.time/30, 2)
    tumor_clin = tumor_clin[!is.na(tumor_clin$ajcc_pathologic_tumor_stage),]
    tumor_clin$stage = tumor_clin$ajcc_pathologic_tumor_stage
    tumor_clin = tumor_clin[grepl("Stage I",tumor_clin$stage),]
    tumor_clin$stage = gsub("[ABCD]$","",tumor_clin$stage)
    tpm_data = tpm_data[, colnames(tpm_data) %in% rownames(tumor_clin)]
    
    # save file for survival analysis
    print(paste0("There are ", ncol(tpm_data), " patients used for analysis!"))
    fastSave::save.pigz(tpm_data, tumor_clin, file = paste0("res/", dir.used, tumor, "_data_for_survival.rda"), n.cores = 10)
    
    return(tumor_clin)
})
all_clin <- do.call(rbind, all_clin)
fastSave::saveRDS.pigz(all_clin, file = paste0("res/", dir.used, "clinical_info_of_5_cancers.rds"))
