genes <- c("Ednrb", "Gfra1", "Ikbkap", "Sox10", "Zeb2", 
           "Fam213a", "Kif1bp", "Phox2b", "Tcf4", "Ubr4", 
           "Ece1", "L1cam", "Nrg1", "Sema3c", "Sema3d", 
           "Acss2", "Ret", "Adamts17", "Eno3", "Sh3pxd2a", 
           "Slc27a4", "Edn3", "Gdnf", "Nrtn", "Gata2", 
           "Rarb", "Nkx2-5", "Cbl")


heatmap_data <- heatmap_data[order.dendrogram(as.dendrogram(hclust(dist(heatmap_data)))),order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data)))))]
rownames(heatmap_data) <- lookupGeneName(dat.filtered, rownames(heatmap_data))

heatmap_annotation_row <- heatmap_annotation_row[,order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_annotation_row)))))]
rownames(heatmap_annotation_row) <- lookupGeneName(dat.filtered, glia_genotype_diff_test_sig_genes)

heatmap_annotation_col <- pData(dat.filtered)[,c("genotype", "age", "cluster")]

heatmap_colors <- list(
  "genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"),
  "cluster" = c("1" = "#66C2A5",
                "2" = "#FC8D62",
                "3" = "#8DA0CB",
                "4" = "#E78AC3",
                "5" = "#A6D854",
                "6" = "#FFD92F",
                "7" = "#E5C494",
                "8" = "#B3B3B3",
                "9" = "#80B1D3"),
  "age" = c("E12.5" = "#4daf4a",
            "E14.5" = "#984ea3"),
  
  "CC_regulation" = c("1" = "#E7298A", "0" = "grey90")
)

#colnames(heatmap_annotation_col) <- c("Cell Type","Genotype")

#pdf("DEG_glia!cluster9_genotype_heatmap_test.pdf", height = 20, width = 20)
pheatmap(mat = heatmap_data,
         scale="none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         drop_levels = T,
         cluster_rows = F,
         cluster_cols = F,
         color = viridis(100),
         annotation_col = heatmap_annotation_col,
         annotation_row = heatmap_annotation_row[,"CC_regulation", drop = F],
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = T
)
dev.off()



################################################################################################

violindat <- merge(pData(dat.filtered)[,c("genotype", "cluster")],
                   #                   t(exprs(dat.filtered)[names(hoxDEG),]),
                   t(exprs(dat.filtered)[names(antposDEG),]),
                   by = 0)

rownames(violindat) <- violindat$Row.names
violindat <- violindat[,-1]

melted_violindat <- melt(violindat, id.vars = c("genotype", "cluster"))

plot_list = list()

#for(i in seq_along(antposDEG)){
for(i in 1:4){
  p = melted_violindat[!melted_violindat$cluster == 9 & melted_violindat$variable == lookupGeneId(dat.filtered, c("Gfra1", "Rarb", "Sox10", "Gata2"))[i],] %>%
    ggplot(., aes(x = genotype, y = value, fill = genotype)) +
    geom_violin(trim = F, scale = "width") +
    stat_summary(fun.y = mean, geom = "point") +
    #geom_jitter() +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none", 
          legend.title = element_text(face = "plain"),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(c("Gfra1", "Rarb", "Sox10", "Gata2")[i])
  plot_list[[i]] = p
}

grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]])








HSCR_genes_subset <- c("Gfra1", "Nrg1", "Phox2b", "Ret", "Sema3c", "Sema3d", "Sox10")
HSCR_genes_subset <- sort(HSCR_genes_subset)

violindat3 <- merge(pData(dat.filtered)[,c("genotype", "cluster")],
                    t(exprs(dat.filtered)[rownames(fData(dat.filtered)[fData(dat.filtered)$gene_short_name %in% HSCR_genes_subset,]),]),
                    by = 0)

rownames(violindat3) <- violindat3$Row.names
violindat3 <- violindat3[,-1]

melted_violindat3 <- melt(violindat3, id.vars = c("genotype", "cluster"))

plot_list = list()

for(i in seq_along(HSCR_genes_subset)){
  p = melted_violindat3[melted_violindat3$variable == lookupGeneId(dat.filtered, HSCR_genes_subset)[i],] %>%
    ggplot(., aes(x = genotype, y = log10(value + 1), fill = genotype)) +
    facet_wrap(~cluster, nrow = 1) +
    geom_violin(trim = T, scale = "width") +
    stat_summary(fun.y = mean, geom = "point", size = 1) +
    scale_color_manual(values = c(brewer.pal(8, "Set2"), brewer.pal(5, "Set3")[5])) +
    scale_fill_manual(values = c("gray70", "darkred")) +
    scale_y_continuous(HSCR_genes_subset[i], expand = c(0,0), breaks = c(0, 0.5, 1)) +
    theme(legend.position = "none", 
          plot.title = element_text(face = "plain", size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    )
  plot_list[[i]] = p
}

pdf("HSCR_genes_subset_violin_log=T_trim=T_lockscale.pdf", height = 11, width = 9, useDingbats = F)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
             plot_list[[5]], plot_list[[6]], plot_list[[7]], ncol = 1
)
dev.off()

pdf("F4A-RetGRN-DEG-UMAPs-horizontal.pdf", height = 4, width = 8, useDingbats = F)
plotUMAP(dat.filtered, markers = HSCR_genes_subset, scaled = T, size = 0.2, nrow = 2) +
  scale_color_viridis(option = "inferno", breaks = seq(0, 2, by = 0.5)) +
  #guides(color = guide_legend(title.position = "top")) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, vjust = -2),
        legend.direction = "horizontal",
        legend.position = c(0.75, 0.2),
        legend.key.height = unit(1, "line"),
        legend.key.width = unit(1.5, "line"))
dev.off()
