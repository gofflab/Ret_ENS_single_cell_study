pdf("UMAP_GFRA1_GFRA2_RET_CFP.pdf", height = 5.5, width = 5.5, useDingbats = F)
plotUMAP(dat.filtered, markers = c("tCFP", "Ret", "Gfra1", "Gfra2"), scaled = T, size = 0.2) +
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

pdf("UMAP_GFRA2_CFP.pdf", height = 3, width = 5.5, useDingbats = F)
plotUMAP(dat.filtered, markers = c("tCFP", "Gfra2"), scaled = T, size = 0.2) +
  scale_color_viridis("log10(CPC + 1)\nscaled to max", option = "inferno", breaks = c(0,1)) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = c(0.85, 0.25),
        legend.key.height = unit(0.7, "line"),
        legend.key.width = unit(0.7, "line"))
dev.off()

pdf("UMAP_GFRA2_CFP_neurons.pdf", height = 2, width = 4, useDingbats = F)
plotUMAP(dat.filtered[,pData(dat.filtered)$cluster %in% c(6:8)], markers = c("tCFP", "Gfra2"), scaled = T, size = 0.4) +
  scale_color_viridis("log10(CPC + 1)\nscaled to max", option = "inferno", breaks = c(0,1)) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = c(0.85, 0.25),
        legend.key.height = unit(0.7, "line"),
        legend.key.width = unit(0.7, "line"))
dev.off()


tmp <- pData(dat.filtered)[,c("UMAP1", "UMAP2")]
tmp$genotype <- ifelse(pData(dat.filtered)$cluster == 8, as.character(pData(dat.filtered)$genotype), "blank")

pdf("UMAP_cluster8.pdf", width = 3, height = 2.5, useDingbats = F)
ggplot(tmp, aes(x = UMAP1, y = UMAP2, color = genotype)) + 
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey90", brewer.pal(3, "Set1")[1:2])) +
  theme(aspect.ratio = 1, 
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10)
  )
dev.off()

tmp$cluster <- ifelse(pData(dat.filtered)$cluster %in% 7:8, pData(dat.filtered)$cluster, 0)
tmp$cluster <- as.factor(tmp$cluster)

pdf("UMAP_cluster7vs8.pdf", width = 3, height = 2.5, useDingbats = F)
ggplot(tmp, aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey95", brewer.pal(8, "Set2")[7:8])) +
  theme(aspect.ratio = 1, 
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10)
  )
dev.off()

tmp <- pData(dat.filtered)[pData(dat.filtered)$cluster %in% 6:8,c("UMAP1", "UMAP2", "genotype", "cluster")]
tmp$cluster8_genotype <- ifelse(tmp$cluster == 8, as.character(tmp$genotype), "blank")

pdf("UMAP_cluster8_neurons.pdf", height = 2, width = 2, useDingbats = F)
ggplot(tmp, aes(x = UMAP1, y = UMAP2, color = cluster8_genotype)) + 
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey95", brewer.pal(3, "Set1")[1:2])) +
  theme(aspect.ratio = 1, 
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10)
  )
dev.off()

tmp$cluster78 <- ifelse(tmp$cluster %in% 7:8, tmp$cluster, 0)
tmp$cluster78 <- as.factor(tmp$cluster78)

pdf("UMAP_cluster7vs8_neurons.pdf", width = 2, height = 2, useDingbats = F)
ggplot(tmp, aes(x = UMAP1, y = UMAP2, color = cluster78)) + 
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey95", brewer.pal(8, "Set2")[7:8])) +
  theme(aspect.ratio = 1, 
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10)
  )
dev.off()

pdf("UMAP_GFRA1_CART_RET_CFP.pdf", height = 5.5, width = 5.5, useDingbats = F)
plotUMAP(dat.filtered, markers = c("tCFP", "Ret", "Gfra1", "Cartpt"), scaled = T, size = 0.2) +
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

pdf("UMAP_F3_IHC_genes_neurons.pdf", height = 4, width = 6, useDingbats = F)
plotUMAP(dat.filtered[,pData(dat.filtered)$cluster %in% 6:8], markers = c("tCFP", "Ret", "Gfra1", "Gfra2", "Cartpt"), scaled = T, size = 0.2) +
  scale_color_viridis("log10(CPC + 1)\nscaled to max", option = "inferno", breaks = c(0,1)) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = c(0.75, 0.25),
        legend.key.height = unit(0.7, "line"),
        legend.key.width = unit(0.7, "line"))
dev.off()

my_color_scale <- scale_color_gradientn("", colors = c("#F1ECB9", "#BEDEE5", "#25C6FE", "#9232B4", "#2C0404"), breaks = c(0, 1), labels = c("min", "max"))
my_theme <-   theme(aspect.ratio = 1,
                    strip.background = element_blank(),
                    strip.text = element_text(size = 20, margin = margin(b = 5)),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    axis.title = element_blank())
#Branch A
pdf("BranchA_markers.pdf", height = 6, width = 7, useDingbats = F)
plotUMAP(dat.filtered, markers = c("Gal", "Vip", "Nos1", "Etv1"), size = 0.5, scale = T) +
  my_color_scale +
  my_theme
dev.off()

#Branch B
pdf("BranchB_markers.pdf", height = 6, width =7, useDingbats = F)
plotUMAP(dat.filtered, markers = c("Bnc2", "Dlx5", "Ndufa4l2", "Mgat4c"), size = 0.5, scale = T) +
  my_color_scale +
  my_theme
dev.off()