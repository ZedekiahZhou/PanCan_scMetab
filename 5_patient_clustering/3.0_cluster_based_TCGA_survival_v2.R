# ==================================================================================================
# Author: Zhou Zhe
# Version: 1.0
# Date: Aug 9, 2022
# ==================================================================================================

library(dplyr)
library(GSVA)
library(survival)
library(survminer)
library(meta)
library(ggsci)
dir.used <- "6_patient_clustering/CluSurv/"
dir.create(file.path("res/", dir.used), recursive = T)
dir.create(file.path("plot/", dir.used), recursive = T)

# I. Load data ------------------------------------------------------------------------------------

dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets
celltypes <- c("Malignant", "Myeloid", "T cells", "B cells", "Endothelial", "Fibroblasts")
cluster_sigs <- lapply(celltypes, function(used_celltype) {
    markers_celltype <- read.csv(file = paste0("res/6_patient_clustering/Cluster_DEG/", used_celltype, "_celltype_DEG.csv"))
    markers_celltype <- subset(markers_celltype, p_val_adj < 0.01)
    markers_celltype <- markers_celltype %>%
        group_by(avg_log2FC > 0) %>%
        slice_max(order_by = abs(avg_log2FC), n = 50)
    
    markers_sam <- read.csv(file = paste0("res/6_patient_clustering/Cluster_DEG/", used_celltype, "_sample_DEG.csv"))
    markers_sam <- subset(markers_sam, p_val_adj < 0.05)
    markers_sam <- markers_sam %>%
        group_by(avg_log2FC > 0) %>%
        slice_max(order_by = abs(avg_log2FC), n = 50)    
    
    res = list(celltype_C1 = markers_celltype$X[markers_celltype$avg_log2FC > 0], 
               celltype_C2 = markers_celltype$X[markers_celltype$avg_log2FC < 0], 
               sam_C1 = markers_sam$X[markers_sam$avg_log2FC > 0], 
               sam_C2 = markers_sam$X[markers_sam$avg_log2FC < 0])
    names(res) = paste0(sub(" ", "-", used_celltype), "_", names(res))
    return(res)
})
cluster_sigs <- do.call(c, cluster_sigs)


# II. GSVA Score ----------------------------------------------------------------------------------
cluster.gsva <- lapply(unique(dataset$tumorType), function(cancer_type) {
    load(paste0("res/6_patient_clustering/SurvPre/", cancer_type, "_data_for_survival.rda"))
    tpm_data <- as.matrix(tpm_data)
    tmp.res <- gsva(tpm_data, cluster_sigs, method = "gsva", kcdf = "Gaussian", mx.diff = T, parallel.sz = 10)
    return(t(tmp.res))
})
names(cluster.gsva) = unique(dataset$tumorType)
saveRDS(cluster.gsva, file = "res/6_patient_clustering/SurvPre/cluster_sig_TCGA_gsva.rds")


# III. TCGA Grouping -----------------------------------------------------------------------------------
cluster.gsva <- readRDS("res/6_patient_clustering/SurvPre/cluster_sig_TCGA_gsva.rds")
group_cutoff = c(0.4, 0.6)
cluster.group <- lapply(cluster.gsva, function(gsva.df) {
    group.df = sapply(as.data.frame(gsva.df), function(x) {
        ifelse(x < quantile(x, group_cutoff[1]), "low", 
               ifelse(x > quantile(x, group_cutoff[2]), "high", NA))
    })
    group.df <- data.frame(group.df)
    tmp_group <- unique(sub("_C[12]", "", colnames(gsva.df)))
    tmp_group.df <- sapply(tmp_group, function(x) {
        tmp_dist = gsva.df[, paste0(x, "_C1")] - gsva.df[, paste0(x, "_C2")]
        ifelse(tmp_dist > 0.05, "C1", 
               ifelse(tmp_dist < -0.05, "C2", NA))
    })
    group.df <- cbind(group.df, tmp_group.df)
})


# IV. Cox analysis --------------------------------------------------------------------------------
all_clin <- readRDS("res/6_patient_clustering/SurvPre/clinical_info_of_5_cancers.rds")
all_clin$cancerType <- ifelse(ifelse(all_clin$cancerType == "COAD" | all_clin$cancerType == "READ", "COAD", 
                                     all_clin$cancerType))
tumors = unique(dataset$tumorType)[1:5]
cluster.survival <- lapply(tumors, function(cancer_type) {
    print(cancer_type)
    tumor_clin <- subset(all_clin, cancerType == cancer_type)
    
    res = lapply(colnames(cluster.group[[cancer_type]]), function(x) {
        tmp_group = cluster.group[[cancer_type]][, x, drop = F]
        tmp_group = tmp_group[rownames(tumor_clin), ]
        dat = tumor_clin
        dat$group = tmp_group
        dat = dat[!is.na(dat$group), c("patient", "gender", "age", "stage", "cancerType", "OS", "OS.time.month", "group")]
        
        if (length(unique(dat$gender))>1) {
            cox = coxph(Surv(OS.time.month,OS) ~ group + gender + age + stage, data=dat)
        } else {
            cox = coxph(Surv(OS.time.month,OS) ~ group + age + stage, data=dat)
        }
        cox.summary = summary(cox)
        nvar = length(unique(dat$group)) - 1
        #
        HR = round(cox.summary$conf.int[1:nvar,1], 2)
        HR.lower = round(cox.summary$conf.int[1:nvar,3], 2)
        HR.upper = round(cox.summary$conf.int[1:nvar,4], 2)
        HR.range = sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
        coef = cox.summary$coefficients[1:nvar,1]
        coef.se = cox.summary$coefficients[1:nvar,3]
        Pval = round(cox.summary$coefficients[1:nvar,5], 4)
        group = gsub("group","",rownames(cox.summary$conf.int)[1:nvar])
        return(data.frame(cancer_type = cancer_type, clu_method = x, group=group, 
                          HR=HR, HR.range=HR.range, coef=coef, coef.se=coef.se, Pval=Pval))
    })
    res = do.call(rbind, res)
    rownames(res) = colnames(cluster.group[[cancer_type]])
    return(res)
})
cluster.survival <- do.call(rbind, cluster.survival)
cluster.survival$padj = p.adjust(cluster.survival$Pval, method = "BH")
write.csv(cluster.survival, file = "res/", dir.used, "cluster.survival.csv", row.names = F)


tmp <- split(cluster.survival, f = cluster.survival$clu_method)
pdf("plot/", dir.used, "cluster_based_TCGA_survival_forest.pdf", width = 14, height = 10)
a = lapply(tmp, function(x) {
    this_method = unique(x$clu_method)
    meta <- metagen(TE = x$coef, seTE = x$coef.se, studlab = x$cancer_type, exclude = x$cancer_type == "PRAD",
                    comb.fixed = F, comb.random = T, prediction = F, sm = "HR")
    
    forest(meta, test.overall.random = T, digits.pval = 4, 
           colgap.forest.left = "5cm", zero.pval = T, xlab = this_method)
    return(meta)
})
dev.off()

b <- reshape2::dcast(cluster.survival, clu_method ~ cancer_type)
rownames(b) <- b$clu_method
b <- b[-1]
b$nDataset <- rowSums(b < 0.05)


pdf("plot/", dir.used, "KM_plot.pdf", width = 5.5, height = 6)
km_plot <- lapply(tumors, function(cancer_type) {
    print(cancer_type)
    tumor_clin <- subset(all_clin, cancerType == cancer_type)
    
    res = lapply(colnames(cluster.group[[cancer_type]]), function(x) {
        tmp_group = cluster.group[[cancer_type]][, x, drop = F]
        tmp_group = tmp_group[rownames(tumor_clin), ]
        dat = tumor_clin
        dat$group = tmp_group
        dat = dat[!is.na(dat$group), c("patient", "gender", "age", "stage", "cancerType", "OS", "OS.time.month", "group")]
        
        fit = survfit(Surv(OS.time.month, OS) ~ group, data = dat)
        p = ggsurvplot(fit, dat, vlegend.labs=unique(dat$group),
                       surv.median.line="none", pval=T, conf.int=F,
                       risk.table=T, risk.table.y.text.col=T,
                       palette=get_palette("jco", length(unique(dat$group))),
                       legend.title="", ggtheme=theme_bw(), 
                       fontsize = 5, 
                       font.x = 16, font.y = 16, font.tickslab = 16, 
                       font.legend = 16, pval.size = 5, 
                       censor.size = 3, size = 0.5) + 
            xlab("Months") + ggtitle(sprintf("data: %s - %s", cancer_type, x))
        print(p)
        return(p)
    })
    names(res) = colnames(cluster.group[[cancer_type]])
    return(res)
})
dev.off()


