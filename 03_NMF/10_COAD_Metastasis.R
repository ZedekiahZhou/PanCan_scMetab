library(dplyr)
library(GSVA)

### sample info
sam_info <- read.delim("data/NMF_valid/CRC_Metastasis/GSE50760/sam_info.txt", header = F)
colnames(sam_info) <- c("GSM", "Description")
sam_info$Group = sub(" .+", "", sam_info$Description)
sam_info$sample = make.names(sub(".+ ", "", sam_info$Description))

### read expression matrix
fpath = dir("data/NMF_valid/CRC_Metastasis/GSE50760/GSE50760_RAW/")
expr <- lapply(fpath, function(x) {
    read.delim(gzfile(paste0("data/NMF_valid/CRC_Metastasis/GSE50760/GSE50760_RAW/", x)), header = TRUE)
})
expr = Reduce(function(x, y) {
    if (all(x$genes == y$genes)) return(cbind(x, y[, 2, drop = F])) else stop("Genes not consistent!")
}, x = expr)
colnames(expr) = sub("_FPKM", "", colnames(expr))
expr <- expr[!duplicated(expr$genes), ]
rownames(expr) <- expr$genes
expr = log2(expr[-1] + 1)
expr = expr[rowSums(expr) > 0, ]

all(sam_info$sample == colnames(expr))

### calculate MP score 
MP_list <- as.list(readRDS("res/02_NMF/MPs/MP_list.RDS")$Malignant)
res.gsva <- gsva(as.matrix(expr), MP_list, method = "gsva", kcdf = "Gaussian", mx.diff = T, parallel.sz = 10)
res.gsva <- data.frame(t(res.gsva), check.names = FALSE)

fun_ttest <- function(x, idx1, idx2) {
    tres = t.test(x[idx1], x[idx2])
    data.frame(t = tres$statistic, p = tres$p.value)
}

res_m2p = plyr::ldply(res.gsva, fun_ttest, 
                      idx1 = which(sam_info$Group == "metastasized"), 
                      idx2 = which(sam_info$Group == "primary"))
res_p2n = plyr::ldply(res.gsva, fun_ttest, 
                      idx1 = which(sam_info$Group == "primary"), 
                      idx2 = which(sam_info$Group == "normal"))

mmp7_gene <- MP_list$`MMP7`[MP_list$`MMP7` %in% rownames(expr)]
mmp7 <- expr[mmp7_gene, ]
pheatmap::pheatmap(mmp7, scale = "row", cluster_cols = FALSE)
