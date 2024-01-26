density_plot <- function(seu.list, 
                         gene = NULL, 
                         metab_path = NULL, 
                         sig_path = NULL, 
                         types = "1D", 
                         cluster_cols = F, 
                         cluster_rows = F
) {
    require(ggpubr)
    require(ggridges)
    require(RColorBrewer)
    require(patchwork)
    plot.list = lapply(1:nrow(dataset), function(i) {
        tumor <- dataset$DataSets[i]
        path2 = path1 = NULL
        pas <- data.frame(row.names = colnames(seu.list[[tumor]]))
        
        if (!is.null(gene)) {
            gene <- intersect(gene, rownames(seu.list[[tumor]]))
            gene_expr <- FetchData(seu.list[[tumor]], vars = gene, slot = "data")
            pas <- cbind(pas, gene_expr)
            if (is.null(path1)) path1 = path2= gene else path2 = gene
        }
        
        if (!is.null(metab_path)) {
            metab_pas <- data.frame(fread(paste0("res/03_SCENIC/aucell/", tumor, "_metab_aucell.csv")),
                                    row.names = 1, check.names = F)
            metab_pas <- metab_pas[colnames(seu.list[[tumor]]), metab_path, drop = FALSE]
            pas <- cbind(pas, metab_pas)
            if (is.null(path1)) path1 = path2= metab_path else path2 = metab_path
        }
        
        if (!is.null(sig_path)) {
            sig_pas <- data.frame(fread(paste0("res/03_SCENIC/aucell/", tumor, "_celltypeSig_aucell.csv")),
                                  row.names = 1, check.names = F)
            sig_pas <- sig_pas[colnames(seu.list[[tumor]]), sig_path, drop = FALSE]
            if (is.null(path1)) path1 = path2= sig_path else path2 = sig_path
            pas = cbind(pas, sig_pas)
        }
        
        
        if (types == "heatmap") {
            res = pheatmap::pheatmap(cor(pas[, path1], pas[, path2]), main = tumor, 
                                      cluster_cols = cluster_cols, cluster_rows = cluster_rows)
        } else if (types == "ridge") {
            df_comb = expand.grid(path1, path2)
            res.list <- apply(df_comb, 1, function(x) {
                pas$group = ifelse(pas[[x[1]]] > 0, paste0(x[1], "+"), paste0(x[1], "-"))
                ggplot(pas, aes(y = group, x = .data[[x[2]]], fill = group)) + 
                    geom_density_ridges() + 
                    scale_fill_manual(values = brewer.pal(9, "Set1")[2:1]) + 
                    theme_bw() + 
                    theme(panel.grid = element_blank())
            })
            res = wrap_plots(res.list, ncol = 4) + plot_annotation(tumor)
        } else {
            # add noise to avoid error of geom_density_2d_filled()
            withr::with_seed(seed = 2021, 
                             noise <- matrix(runif(nrow(pas)*ncol(pas), min = 0, max = 1e-10), ncol = ncol(pas)))
            pas <- pas + noise
            if (types == "1D") {
                res.list <- lapply(colnames(pas), function(x) {
                    ggplot(pas, aes(.data[[x]])) + 
                        geom_density(color = "black", fill = "grey") + 
                        theme_bw() + 
                        theme(legend.position = "none", panel.grid = element_blank())
                })
            } else {
                if (types == "Self-2D") {
                    df_comb = t(combn(path1, m = 2))
                } else if (types == "2D") {
                    df_comb = expand.grid(path1, path2)
                } else stop('"types" should be one of c("heatmap", "1D", "Self-2D", "2D")')
                
                res.list <- apply(df_comb, 1, function(x) {
                    ggplot(pas, aes(.data[[x[1]]], .data[[x[2]]])) + 
                        geom_density_2d_filled(contour_var = "ndensity") + 
                        scale_fill_manual(values = c("white", brewer.pal(9, "YlOrRd"))) + 
                        stat_cor(label.y.npc = 1) +
                        # stat_regline_equation(label.y.npc = 0.95) + 
                        theme_bw() + 
                        theme(legend.position = "none", panel.grid = element_blank())
                })
            }
            res = wrap_plots(res.list, ncol = 4) + plot_annotation(tumor)
        }
        return(res)
    })
    return(plot.list)
}