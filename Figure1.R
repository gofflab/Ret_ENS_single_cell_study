#1A
my_collapse_labels_lines <- function (labels) 
{
  out <- do.call("Map", c(list(paste, sep = " "), labels))
  list(unname(unlist(out)))
}

my_label_value <- function (labels, multi_line = TRUE) 
{
  labels <- lapply(labels, as.character)
  if (multi_line) {
    labels
  }
  else {
    my_collapse_labels_lines(labels)
  }
}

variable_labeller <- function(labels){
my_label_value(labels, multi_line = F)
  }


p1 <- ggplot(pData(dat.filtered), aes(x = genotype)) + 
  geom_bar(aes(fill = genotype)) + 
  facet_wrap(~age + factor(pData(dat.filtered)$sex, labels = c("F", "M")), nrow = 1, labeller = variable_labeller) + 
  scale_y_continuous("# Cells", expand = expand_scale(mult = c(0, .05)), limits = c(0, 186)) +
  scale_fill_brewer("Genotype", palette="Set1") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.4) +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(b = 5)))

pdf("Cells_passing_filters.pdf", width = 3.5, height = 3)
p1
dev.off()

#1B
p2 <- plotUMAP(dat.filtered, color = "genotype", size = 0.1) +
  guides(color = guide_legend(title = "Genotype", override.aes = list(size = 3))) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.5, 0.2),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.key.height = unit(0.8, "line"),
        legend.key.width = unit(0.5, "line")
        )

p2

pdf("UMAP_genotype.pdf", width = 3, height = 2.5, useDingbats = F)
p2
dev.off()

#1C
p3 <- plotUMAP(dat.filtered, color = "cluster", size = 0.1) +
  scale_color_manual(values = c(brewer.pal(8, "Set2"), brewer.pal(5, "Set3")[5])) +
  guides(color = guide_legend(title = "Cluster", ncol = 3, override.aes = list(size = 2.5))) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.4, 0.2),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.key.height = unit(0.8, "line"),
        legend.key.width = unit(0.5, "line"),
        legend.title = element_text(margin = margin(l = 2, b = -4))
  )


p3

pdf("UMAP_cluster.pdf", width = 3, height = 2.5, useDingbats = F)
p3
dev.off()

#1C
p4 <- plotUMAP(dat.filtered, color = "celltype", size = 0.1) +
  scale_color_manual(values = cbPalette[1:2]) +
  guides(color = guide_legend(title = "Cell Type", override.aes = list(size = 2.5))) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.4, 0.2),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.key.height = unit(0.8, "line"),
        legend.key.width = unit(0.5, "line"),
        legend.title = element_text(margin = margin(l = 2, b = -4))
  )


p4

pdf("UMAP_celltype.pdf", width = 3, height = 2.5, useDingbats = F)
p4
dev.off()


# Supp
p5 <- plotUMAP(dat.filtered, color = "age", size = 0.1) +
  scale_color_manual(values = c("#4DAF4A", "#984EA3")) +
  guides(color = guide_legend(title = "Age", override.aes = list(size = 2.5))) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.4, 0.2),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.key.height = unit(0.8, "line"),
        legend.key.width = unit(0.5, "line"),
        legend.title = element_text(margin = margin(l = 2, b = -4))
  )


p5

pdf("UMAP_age.pdf", width = 3, height = 2.5, useDingbats = F)
p5
dev.off()


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
-ggplot(tmp, aes(x = UMAP1, y = UMAP2, color = genotype)) + 
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


pdf("pattern42.pdf", width = 3, height = 2.5, useDingbats = F)
ggplot(tmp2[tmp2$pattern == "Pattern_42" & tmp2$subset == "HD",]) +
  geom_point(aes(x=x, 
                 y=y,
                 color = weight),
             size = 0.1) +
  scale_color_viridis(option ="inferno", guide = F) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 8, face = "plain", hjust = 0.01, margin = margin(b = -16)),
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.title = element_text(size = 8),
        axis.title.x = element_text(margin = margin(t = -1)),
        axis.title.y = element_text(margin = margin(r = -1)),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  ) 
dev.off()


pdf("MPM.pdf", width = 3, height = 2.5, useDingbats = F)
ggplot(CC_patterns$projection) +
  geom_point(aes(x=UMAP1, y=UMAP2, 
                 color = max_pattern,
                 alpha = max_value), size = 0.2,
             shape = 16) +
  scale_color_brewer(palette = "Dark2") +
  scale_alpha(range = c(0,1)) +
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none"
  )

dev.off()
