pdf("SF2A-UMAP_parentcelltype_prediction.pdf", width = 4, height = 3, useDingbats = F)

plotUMAP(dat.filtered, color = "parentCelltype", size = 0.1) +
  scale_color_manual(values = cbPalette[1:2]) +
  guides(color = guide_legend(title = "Cell Type", override.aes = list(size = 2))) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.4, 0.2),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
  )

dev.off()


pdf("SF2B-UMAP_cluster_pred.pdf", width = 6, height = 3, useDingbats = F)

plotUMAP(dat.filtered, color = "celltype_pred", size = 0.1) +
  scale_color_viridis(discrete = T, option = "H", drop = F) +
  guides(color = guide_legend(title = "Predicted Cluster", override.aes = list(size = 2), ncol = 2)) +
  theme(aspect.ratio = 1, 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
  )

dev.off()


pdf("SF2C-Linnarsson_cluster.pdf", width = 8, height = 3, useDingbats = F)

ggplot(linnarsson_dat, aes(x = X, y = Y)) +
  geom_point(aes(color = cluster_name), size = 0.1) +
  facet_wrap(~parentCellType) +
  scale_color_viridis(discrete = T, option = "H", drop = F) +
  guides(color = guide_legend(title = "Cluster", override.aes = list(size = 2), ncol = 2)) +
  scale_x_continuous("Reduced Dimension 1") +
  scale_y_continuous("Reduced Dimension 2") +
  theme(aspect.ratio = 1, 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(b = 3)),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
  )

dev.off()

pdf("SF2-Linnarsson_description.pdf", width = 9, height = 3, useDingbats = F)

ggplot(linnarsson_dat, aes(x = X, y = Y)) +
  geom_point(aes(color = factor(linnarsson_dat$description, 
                                levels = c("Enteric glia, proliferating", "Enteric glia", "Enteric mesothelial fibroblasts", 
                                           "Cholinergic enteric neurons", "Cholinergic enteric neurons, VGLUT2", "Nitrergic enteric neurons"))),
             size = 0.1) +
  facet_wrap(~parentCellType) +
#  scale_color_viridis(discrete = T, option = "C", drop = F) +
  scale_color_brewer(palette = "Set1") +
  guides(color = guide_legend(title = "Description", override.aes = list(size = 2), ncol = 1)) +
  scale_x_continuous("Reduced Dimension 1") +
  scale_y_continuous("Reduced Dimension 2") +
  theme(aspect.ratio = 1, 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(b = 3)),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
  )

dev.off()


pdf("SF2D-cluster_proportion_celltype.pdf", height = 3, width = 3, useDingbats = F)
ggplot(pData(dat.filtered)) + 
  geom_bar(aes(x = cluster, fill = parentCelltype), position = "fill") + 
  geom_hline(aes(yintercept = 0.1), linetype = "dashed", col = "grey20", size = 1) +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed", col = "grey20", size = 1) +
  scale_fill_manual(values = cbPalette[1:2], guide = F) +
  scale_y_continuous("Proportion", expand = c(0,0), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)) +
  scale_x_discrete("Cluster") +
  theme(plot.margin = margin(6, 0, -10, 0),
        plot.background = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 2),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = 4),
        axis.title.y = element_text(vjust = 0)) +
  monocle:::monocle_theme_opts()
dev.off()

