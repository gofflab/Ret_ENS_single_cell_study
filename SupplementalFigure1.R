callouts <- list()
callouts[[1]] <- c("Myl9", "Plp1", "Fabp7", "Ednrb")
callouts[[2]] <- c("Hist1h2ah", "Hist1h2ai", "Hist1h2ap", "Hist1h1b", "Hist1h2ag", "Hist1h2ao")
callouts[[3]] <- c("Aurka", "Aurkb", "Top2a", "Cdk1", "Ube2c")
callouts[[4]] <- c("Nfix", "Neurod4")
callouts[[5]] <- c("Pcp2", "H2-Ab1", "Hes5")
callouts[[6]] <- c("Hes6", "Btg2")
callouts[[7]] <- c("Nos1", "Ache", "Vip", "Etv1", "Gal", "Cartpt")
callouts[[8]] <- c("Ndufa4l2", "Mgat4c", "Bnc2", "Dlx5", "Npy", "Gfra2", "Slc18a3", "Snap25") 
callouts[[9]] <- c("Hoxc10", "Hoxa10", "Hoxc6", "Hoxa7", "Hoxc9", "Hoxa9", "Gata2", "Gfra3", "Hoxb7", "Phox2a")

plot_specificity <- function(x){
  #  ggplot(melted_df[melted_df$cluster == x,], aes(x = log10(mean_expr + 1), 
  df <- melted_df[melted_df$cluster == x & melted_df$gene_short_name %in% gene_specificity[[x]],]
  ggplot(df, aes(x = log10(mean_expr + 1), y = specificity)) +
    ggtitle(paste("Cluster", x)) +
    geom_point(aes(color = ifelse(ifelse(df$gene_short_name %in% callouts[[x]],
                                         df$gene_short_name, 
                                         "") == "", "no", "yes"),
                   alpha = ifelse(ifelse(df$gene_short_name %in% callouts[[x]],
                                         df$gene_short_name, 
                                         "") == "", "no", "yes")),
               size = 0.1) +
    scale_color_manual(values = c("no" = "black", "yes" = "red")) +
    scale_alpha_discrete(range = c(0.4, 1)) +
    scale_x_continuous("log10(mean expression + 1)", expand = c(0,0)) +
    scale_y_continuous("Specificity", expand = c(0,0), limits = c(0,1)) +
    geom_text(aes(label = ifelse(df$gene_short_name %in% callouts[[x]],
                                 df$gene_short_name, 
                                 "")),
              size = 4, 
              hjust = 0, 
              vjust = 0, 
              nudge_x = 0.02, 
              nudge_y = -0.004, 
              color = "black") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), size = 0.4) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, hjust = 0.005, face = "plain"),
          plot.margin = margin(0, 2, 2, 2),
          axis.title = element_text(size= 10),
          axis.title.x = element_text(margin = margin(t = 0)),
          axis.title.y = element_text(margin = margin(r = 0)),
          axis.text = element_text(size = 10),
          aspect.ratio = 1
    )
}


p1 <- plot_specificity(1)
p2 <- plot_specificity(2)
p3 <- plot_specificity(3)
p4 <- plot_specificity(4)
p5 <- plot_specificity(5)
p6 <- plot_specificity(6)
p7 <- plot_specificity(7)
p8 <- plot_specificity(8)
p9 <- plot_specificity(9)

pdf("SF1A-cluster_specificity_callouts.pdf", height = 12, width = 12, useDingbats = F)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
dev.off()






pdf("SF1B-cluster_markers.pdf", height = 5.5, width = 5.5, useDingbats = F)
plotUMAP(dat.filtered, 
         markers = c("Plp1", "Myl9", "Aurka", "Aurkb", "Cdk1", "Hes6", "Btg2", 
                     "Etv1", "Vip", "Gal", "Cartpt", "Slc18a3", "Npy", "Snap25", 
                     "Gfra2", "Dbh", "Th", "Ddc"), scaled = T, size = 0.2) +
  scale_color_viridis("log10(CPC + 1)\nscaled to max", option = "inferno", breaks = c(0,1)) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = c(0.85, 0.1),
        legend.key.height = unit(0.7, "line"),
        legend.key.width = unit(0.7, "line"))
dev.off()

markers = c("Plp1", "Myl9", "Aurka", "Aurkb", "Cdk1", "Hes6", "Btg2", 
            "Etv1", "Vip", "Gal", "Cartpt", "Slc18a3", "Npy", "Snap25", 
            "Gfra2", "Dbh", "Th", "Ddc")
genes <- exprs(dat.filtered)[rownames(fData(dat.filtered)) %in% lookupGeneId(dat.filtered, markers), , drop = F]
genes <- log10(genes + 1)
geneMax <- rowMax(genes)
genes <- genes/geneMax
genes <- t(genes)
genes <- melt(genes)
colnames(genes) <- c("cell_id", "gene_id", "value")
genes <- merge(genes, fData(dat.filtered)[, c("gene_id", "gene_short_name")], by.x = "gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
tmp <- merge(pData(dat.filtered)[,c("UMAP1", "UMAP2")], genes, by.x = 0, by.y = "cell_id", sort = FALSE)
tmp$gene_short_name <- factor(tmp$gene_short_name, levels = markers)

pdf("SF1B-cluster_markers.pdf", width = 10, height = 10, useDingbats = F)
ggplot(tmp, aes(x = UMAP1, y = UMAP2)) + 
  geom_point(aes(color = value), size = 0.1) +
  scale_color_viridis("log10(CPC + 1)\nscaled to max", option = "inferno", breaks = c(0,1)) +
  facet_wrap('gene_short_name') +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = c(0.8, 0.1),
        legend.key.height = unit(0.8, "line"),
        legend.key.width = unit(0.8, "line"))
dev.off()



markers = c("Ngfr", "Tubb3")
genes <- exprs(dat.filtered)[rownames(fData(dat.filtered)) %in% lookupGeneId(dat.filtered, markers), , drop = F]
genes <- log10(genes + 1)
geneMax <- rowMax(genes)
genes <- genes/geneMax
genes <- t(genes)
genes <- melt(genes)
colnames(genes) <- c("cell_id", "gene_id", "value")
genes <- merge(genes, fData(dat.filtered)[, c("gene_id", "gene_short_name")], by.x = "gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
tmp <- merge(pData(dat.filtered)[,c("UMAP1", "UMAP2")], genes, by.x = 0, by.y = "cell_id", sort = FALSE)
tmp$gene_short_name <- factor(tmp$gene_short_name, levels = markers)

pdf("SF1C-canonical_markers.pdf", width = 5, height = 4, useDingbats = F)
ggplot(tmp, aes(x = UMAP1, y = UMAP2)) + 
  geom_point(aes(color = value), size = 0.2) +
  scale_color_viridis("log10(CPC + 1)\nscaled to max", option = "inferno", breaks = c(0,1)) +
  facet_wrap('gene_short_name') +
  ggtitle("Canonical Markers") +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
dev.off()
