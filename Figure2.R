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

plotUMAP(dat.filtered, markers = c("Zfp804a"), size = 0.5, scale = T) +
  my_color_scale +
  my_theme


plotUMAP(dat.filtered, markers = c("Hes1", "Sox10", "Sox5", "Sox8"), size = 0.5, scale = T) +
  my_color_scale
