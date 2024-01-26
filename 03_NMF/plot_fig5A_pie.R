rm(list=ls())

library(ggplot2)
library(RColorBrewer)

sam_info_list <- readRDS("data/rds/NMF/cells_info_per_sam.RDS")
sam_info_list <- sam_info_list[c(-1, -6)]

MP_list <- readRDS("res/02_NMF/MPs/MP_list.RDS")
MP_list <- MP_list[-1]

df = data.frame(celltype = factor(names(sam_info_list), levels = names(sam_info_list)), 
                n_samples = sapply(sam_info_list, function(x) nrow(x$sam_info)), 
                n_cells = sapply(sam_info_list, function(x) length(unlist(x$sam_cells))), 
                n_MPs = sapply(MP_list, function(x) ncol(x)))

mycolor = c(scMetab::celltype_color[c(3, 5)], brewer.pal(8, "Dark2")[3], 
            brewer.pal(12, "Paired")[10], scMetab::celltype_color[c(6, 7)])
names(mycolor) = df$celltype
ggplot(df, aes(x = "", y = n_cells, fill = celltype)) + 
    geom_bar(stat="identity", width=1, color="white") +
    scale_fill_manual(values = mycolor) +
    coord_polar("y", start=0) +
    theme_void()
ggsave("plot/02_NMF/Fig5A_TME_pie.pdf", width = 7, height = 7)
