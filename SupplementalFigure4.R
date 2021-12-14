tmp <- merge(pData(dat.filtered)[,c("UMAP1", "UMAP2")], np25.gapsResult@sampleFactors, by = 0)
plot_list2 = list()
for (i in 1:25) {
  if(i == 1){
    mytheme <- theme(plot.title = element_text(hjust = 0.005, vjust = -5, size = 10, face = "plain", family = "Helvetica"),
                     aspect.ratio = 1,
                     axis.title = element_text(size = 8),
                     axis.title.x = element_text(margin = margin(t = -2)),
                     axis.title.y = element_text(margin = margin(r = -1)),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     legend.position = c(0.75, 0.25), 
                     legend.title = element_text(size = 6),
                     legend.text = element_text(size = 6),
                     legend.key.size = unit(0.25, "cm"))
  }else{
    mytheme <- theme(plot.title = element_text(hjust = 0.005, vjust = -5, size = 10, face = "plain", family = "Helvetica"),
                     aspect.ratio = 1,
                     axis.title = element_text(size = 8),
                     axis.title.x = element_text(margin = margin(t = -2)),
                     axis.title.y = element_text(margin = margin(r = -1)),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     legend.position = "none")
  }
  p = ggplot(tmp) +
    geom_point(aes_string(x="UMAP1", y="UMAP2", color = paste0("Pattern_", i)), size = 0.1) +
    scale_color_viridis("Pattern\nweight", option = "inferno", breaks = c(0, 0.5, 1), limits = c(0,1)) +
    ggtitle(paste("Pattern", i)) +
    mytheme
  plot_list2[[i]] = p
}

pdf("np25_multipage.pdf", height = 3, width = 3, onefile = T, useDingbats = F)
for (i in seq_along(plot_list2)) {
  print(plot_list2[[i]])
}
dev.off()


pdf("np25_singlepage.pdf", width = 12, height = 12, useDingbats = F)
do.call("grid.arrange", c(plot_list2, ncol = 5))
dev.off()
