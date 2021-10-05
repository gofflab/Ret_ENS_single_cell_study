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
