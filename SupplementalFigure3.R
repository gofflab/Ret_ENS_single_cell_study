pdf("SF3_UMAP_neuron_branch_markers.pdf", height = 4, width = 7, useDingbats = F)
plotUMAP(dat.filtered[,pData(dat.filtered)$cluster %in% 6:8], 
         markers = c("Etv1", "Nos1", "Vip", "Gal", "Bnc2", "Dlx5", "Mgat4c", "Ndufa4l2"), 
         scaled = T, size = 0.2, nrow = 2) +
  scale_color_viridis("log10(CPC + 1)\nscaled to max", option = "inferno", breaks = c(0,1)) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
#        legend.position = c(0.75, 0.25),
        legend.key.height = unit(0.7, "line"),
        legend.key.width = unit(0.7, "line"))
dev.off()
