# ==================================================================================================
# Author: Zhou Zhe
# Program: TCGA survival analysis for metabolic pathway and TF
# Version: 1.0
# Date: Aug 18, 2022
# ==================================================================================================
rm(list=ls())
library(SummarizedExperiment)
library(data.table)
library(scMetab)
library(survival)
library(survminer)
library(GSEABase)
library(parallel)
library(GSVA)
library(meta)


dataset <- as.data.frame(readxl::read_excel("../Datasets.xlsx"))
rownames(dataset) <- dataset$DataSets

# I. GSVA Score ----------------------------------------------------------------------------------
path.gsva <- lapply(dataset$DataSets[1:9], function(tumor) {
    load(paste0("res/survival/", dataset[tumor, ]$tumorType, "_data_for_survival.rda"))
    tpm_data <- as.matrix(tpm_data)
    
    regulons <- geneIds(getGmt(paste0("res/5_SCENIC/regulon/", tumor, "_reg_metab.gmt")))
    TF_expr <- t(tpm_data[names(regulons), ])
    names(regulons) <- paste0(names(regulons), "_reg")
    metab <- geneIds(getGmt(paste0("res/5_SCENIC/dataset_gmt/", tumor, "_gaudeMetab.gmt")))
    
    gls <- c(metab, regulons)
    
    res.gsva <- gsva(tpm_data, gls, method = "gsva", kcdf = "Gaussian", mx.diff = T, parallel.sz = 10)
    res.gsva <- t(res.gsva)
    
    tmp.res <- cbind(res.gsva, TF_expr[rownames(res.gsva), ])
    return(tmp.res)
})
names(path.gsva) = dataset$DataSets[1:9]
fastSave::saveRDS.pigz(path.gsva, file = "res/6_patient_clustering/SurvPre/path_TF_TCGA_gsva.rds", n.cores = 10)


# II. TCGA grouping -------------------------------------------------------------------------------
path.gsva <- readRDS("res/6_patient_clustering/SurvPre/path_TF_TCGA_gsva.rds")
group_cutoff = c(0.4, 0.6)
path.group <- lapply(path.gsva, function(gsva.df) {
    group.df = sapply(data.frame(gsva.df, check.names = F), function(x) {
        ifelse(x < quantile(x, group_cutoff[1]), "low", 
               ifelse(x > quantile(x, group_cutoff[2]), "high", NA))
    })
    group.df <- data.frame(group.df, check.names = F, row.names = rownames(gsva.df))
})


# III. Cox analysis -------------------------------------------------------------------------------
all_clin <- readRDS("res/6_patient_clustering/SurvPre/clinical_info_of_5_cancers.rds")
all_clin$cancerType <- ifelse(ifelse(all_clin$cancerType == "COAD" | all_clin$cancerType == "READ", "COAD", 
                                     all_clin$cancerType))

path.survival <- lapply(dataset$DataSets[1:9], function(tumor) {
    cancer_type <- dataset[tumor, ]$tumorType
    print(paste(tumor, "of", cancer_type))
    tumor_clin <- subset(all_clin, cancerType == cancer_type)
    
    res = lapply(colnames(path.group[[tumor]]), function(x) {
        # print(x)
        tmp_group = path.group[[tumor]][, x, drop = F]
        tmp_group = tmp_group[rownames(tumor_clin), ]
        dat = tumor_clin
        dat$group = tmp_group
        dat = dat[!is.na(dat$group), c("patient", "gender", "age", "stage", "cancerType", "OS", "OS.time.month", "group")]
        
        # some features can not group sample into 2 group
        if (length(unique(dat$group)) <= 1) {
            return(data.frame(cancer_type = cancer_type, dataset = tumor, 
                              clu_method = x, group=NA, 
                              HR=NA, HR.range=NA, coef=NA, coef.se=NA, Pval=NA))
        }
        
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
        return(data.frame(cancer_type = cancer_type, dataset = tumor, 
                          clu_method = x, group=group, 
                          HR=HR, HR.range=HR.range, coef=coef, coef.se=coef.se, Pval=Pval))
    })
    res = do.call(rbind, res)
    rownames(res) = colnames(path.group[[tumor]])
    res$padj <- p.adjust(res$Pval, method = "BH")
    return(res)
})
path.survival <- do.call(rbind, path.survival)
write.csv(path.survival, file = "res/6_patient_clustering/PathSurv/path.survival.csv", row.names = F)


### Plot pathways ---------------------------------------------------------------------------------
path.survival = read.csv("res/6_patient_clustering/PathSurv/path.survival.csv")
plot_data <- subset(path.survival, dataset %in% c("BRCA", "LUAD", "COAD", "PAAD", "STAD") &
                        clu_method %in% gaude_trans$original)
plot_data$clu_method <- gaude_trans[plot_data$clu_method, ]$updated
idx <- plot_data$Pval > 0.05
plot_data$Survival <- ifelse(plot_data$Pval > 0.05, "No_effect", 
                             ifelse(plot_data$HR < 1, "Worse", "Better"))
used_path <- as.data.frame(sort(table(plot_data$clu_method[plot_data$Pval < 0.05]), decreasing = T))
used_path <- subset(used_path, used_path$Freq > 1)
plot_data <- subset(plot_data, clu_method %in% used_path$Var1)
plot_data$clu_method <- factor(plot_data$clu_method, levels = used_path$Var1)

# tmpfile <- paste0("plot/", dir.used, "ALL_pooled_dotplot.pdf")
tmpcolor <- brewer.pal(9, "Set1")[c(1,2,11)]
pdf("plot/6_patient_clustering/PathSurv/Pathway_survival_dotplot.pdf", width = 6, height = 4)
# pdf(tmpfile, width = 12, height = length(levels(TvN_metab_markers$pathway))/9)
print(ggplot(plot_data, aes(x = clu_method, y = cancer_type, color = Survival)) + 
          geom_point(aes(size = -log10(Pval))) + 
          scale_color_manual(values = c(Worse = tmpcolor[1], Better = tmpcolor[2], No_effect = tmpcolor[3])) + 
          scale_size(limits = c(0, 4), range = c(0, 4)) + 
          theme_classic() + 
          xlab("Pathway") + ylab("Cancer type") + 
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black"), 
                legend.position = "top", 
                plot.margin = margin(0.5, 2, 0.5, 2, unit = "cm")))
dev.off()

pdf("plot/6_patient_clustering/PathSurv/Pathway_survival_dotplot_vertical.pdf", width = 4.6, height = 4.5)
# pdf(tmpfile, width = 12, height = length(levels(TvN_metab_markers$pathway))/9)
print(ggplot(plot_data, aes(y = clu_method, x = cancer_type, color = Survival)) + 
          geom_point(aes(size = -log10(Pval))) + 
          scale_color_manual(values = c(Worse = tmpcolor[1], Better = tmpcolor[2], No_effect = tmpcolor[3])) + 
          scale_size(limits = c(0, 4), range = c(0, 4)) + 
          theme_classic() + 
          ylab("Pathway") + xlab("Cancer type") + 
          theme(axis.text=tmp.text, text = tmp.text, 
                axis.text.x = element_text(angle = 45, hjust = 1), 
                axis.line=element_line(size = .3, colour="black"), 
                legend.position = "top", 
                plot.margin = margin(0.5, 2, 0.5, 2, unit = "cm")))
dev.off()



## IV. Meta analysis -----------------------------------------------------------------------------
# remove data from validate set
path.survival <- read.csv("res/6_patient_clustering/PathSurv/path.survival.csv")
path.survival <- subset(path.survival, dataset %in% c("BRCA", "COAD", "LUAD", "PAAD", "STAD"))

tmp <- split(path.survival, f = path.survival$clu_method)
nDataset <- sapply(tmp, nrow)
tmp <- tmp[nDataset == 5]
pdf("plot/6_patient_clustering/PathSurv/Path_TCGA_survival_forest_Bigger_LungColon.pdf", width = 14, height = 10)
a = lapply(tmp, function(x) {
    this_method = unique(x$clu_method)
    meta <- metagen(TE = x$coef, seTE = x$coef.se, studlab = x$cancer_type, # exclude = x$cancer_type == "PRAD",
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

# too slow, re-run in caution!
pdf("plot/6_patient_clustering/PathSurv/Path_KM_plot_Bigger_LungColon.pdf", width = 5.5, height = 6)
km_plot <- lapply(dataset$DataSets[c(1,4,7:9)], function(cancer_type) {
    print(paste(cancer_type))
    tumor_clin <- subset(all_clin, cancerType == cancer_type)
    
    res = lapply(colnames(path.group[[cancer_type]]), function(x) {
        tmp_group = path.group[[cancer_type]][, x, drop = F]
        tmp_group = tmp_group[rownames(tumor_clin), ]
        dat = tumor_clin
        dat$group = tmp_group
        dat = dat[!is.na(dat$group), c("patient", "gender", "age", "stage", "cancerType", "OS", "OS.time.month", "group")]
        
        # some features can not group sample into 2 group
        if (length(unique(dat$group)) <= 1) {
            return(NULL)
        }
        
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
    names(res) = colnames(path.group[[cancer_type]])
    return(res)
})
dev.off()

