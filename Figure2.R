#Figure 2A
heatmap_data <- log10(exprs(dat.filtered)[genotype_diff_test_sig_genes,] + 1) # log transform
heatmap_data <- t(apply(heatmap_data, 1, function(x){ # Normalize gene expr by dividing by max
  x/max(x)
}))

rownames(heatmap_data) <- lookupGeneName(dat.filtered, rownames(heatmap_data))

heatmap_data <- heatmap_data[order.dendrogram(as.dendrogram(hclust(dist(heatmap_data)))),order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data)))))]

heatmap_annotation <- pData(dat.filtered)[, "genotype", drop = F]
colnames(heatmap_annotation) <- "Genotype"

heatmap_colors <- list(
  "Genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"))

pdf("DEG_all_cells_genotype_heatmap.pdf", height = 20, width = 20)
pheatmap(mat = heatmap_data,
         scale="none",
         show_rownames = TRUE,
         show_colnames = FALSE,
         drop_levels = FALSE,
         cluster_rows = F,
         cluster_cols = F,
         color = viridis(100),
         annotation_col = heatmap_annotation,
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         fontsize_row = 5
)
dev.off()

#Figure 2C

hoxDEG <- c("Hoxa2", "Hoxb2", "Hoxd3", "Hoxa4", "Hoxa5", "Hoxb5", "Hoxb6")
names(hoxDEG) <- lookupGeneId(dat.filtered, hoxDEG)

violindat <- merge(pData(dat.filtered)[,c("genotype", "cluster")],
                   t(exprs(dat.filtered)[names(hoxDEG),]),
                   by = 0)
rownames(violindat) <- violindat$Row.names
violindat <- violindat[,-1]

melted_violindat <- melt(violindat, id.vars = c("genotype", "cluster"))

plot_list = list()

for(i in seq_along(hoxDEG)){
  p = melted_violindat[melted_violindat$variable == names(hoxDEG)[i],] %>%
    ggplot(., aes(x = genotype, y = log10(value + 1), color = genotype, fill = genotype)) +
    geom_violin() +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none") +
    ggtitle(hoxDEG[i])
  plot_list[[i]] = p
}

pdf("F2C-hoxDEG-by-genotype.pdf", width = 5.5, height = 10)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
             plot_list[[4]], plot_list[[5]], plot_list[[6]], 
             plot_list[[7]], ncol = 2)
dev.off()
