linnarsson_gapsresult <- readRDS(file = "/Users/liz/Documents/Hopkins/Goff Lab/Chakravarti/Linnarsson_enteric_data/Linnarsson_np70.rds")
linnarsson_dat <- pData(readRDS(file = "../Linnarsson_enteric_data/dat.RDS"))

mat <- log10(exprs(dat.filtered[expressed_genes,]) + 1)
rownames(mat) <- gsub("\\..*", "", rownames(mat))

projection70 <- projectR(mat, linnarsson_gapsresult@featureLoadings, AnnotionObj = NA, IDcol=NA,full=TRUE)

# Plot Linnarsson patterns with HSCR projections 
tmp <- melt(linnarsson_gapsresult@sampleFactors)
tmp1 <- merge(tmp, linnarsson_dat[,c("X", "Y","parentCellType")], by.x = "Var1", by.y = 0)
names(tmp1) <- c("cell_ID", "pattern", "value", "x", "y", "subset")
tmp1$significant <- TRUE

tmp2 <- apply(t(projection70$projection), 2, function(x){
  x=x-min(x)
  x=x/max(x)
})
tmp2 <- melt(tmp2)
#tmp2 <- melt(pData(dat.filtered)[, c("cell_id", paste0("Linnarsson_Pattern_", c(1:70)))])
tmp3 <- merge(tmp2, pData(dat.filtered)[,c("UMAP1", "UMAP2", "project")], by.x = "Var1", by.y = 0)

names(tmp3) <- c("cell_ID", "pattern", "value", "x", "y", "subset")
tmp4 <- melt(t(projection70$pval) < 0.05/nrow(pData(dat.filtered)))
tmp5 <- merge(tmp3, tmp4, by.x = c("cell_ID", "pattern"), by.y = c("Var1", "Var2"), all = TRUE)
names(tmp5) <- c("cell_ID", "pattern", "value", "x", "y", "subset", "significant")

tmp6 <- rbind(tmp1, tmp5)
tmp6$subset <- as.factor(tmp6$subset)
tmp6$subset <- factor(tmp6$subset, levels(tmp6$subset)[c(1, 3, 2)])

pdf("HSCR_projected_into_Linnarsson_2.pdf", height = 40, width = 40)
ggplot(tmp6[tmp6$pattern %in% paste0("Pattern_", c(1:70)),]) +
  geom_point(aes(x=x, y=y, 
                 color = ifelse(significant, value, NA),
                 alpha = ifelse(significant, 1, 0.1)), size = 0.1) +
  facet_wrap(pattern~subset, scales = "free", ncol = 12) +
  scale_color_viridis("Pattern weight", option = "inferno") +
  scale_alpha(guide = "none") +
  scale_x_continuous("Reduced dimension 1") +
  scale_y_continuous("Reduced dimension 2") +
  theme(strip.background = element_blank(), strip.text = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", legend.key.width = unit(60, "pt"))
dev.off()

plot_list = list()
for (i in 1:70) {
    if(i == 1){
    mytheme <- theme(plot.title = element_text(hjust = 0.005, vjust = -5, size = 10, face = "plain", family = "Helvetica"),
          aspect.ratio = 1,
          strip.background = element_blank(), 
          strip.text = element_blank(), 
          axis.title = element_text(size = 8),
          axis.title.x = element_text(margin = margin(t = -2)),
          axis.title.y = element_text(margin = margin(r = -1)),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = c(0.94, 0.25), 
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.25, "cm"))
    }else{
      mytheme <- theme(plot.title = element_text(hjust = 0.005, vjust = -5, size = 10, face = "plain", family = "Helvetica"),
                       aspect.ratio = 1,
                       strip.background = element_blank(), 
                       strip.text = element_blank(), 
                       axis.title = element_text(size = 8),
                       axis.title.x = element_text(margin = margin(t = -2)),
                       axis.title.y = element_text(margin = margin(r = -1)),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       legend.position = "none")
    }
  p = ggplot(tmp6[tmp6$pattern == paste0("Pattern_", i),]) +
    geom_point(aes(x=x, y=y, color = value), size = 0.1) +
    facet_wrap(pattern~subset, scales = "free", ncol = 3) +
    scale_color_viridis("Pattern\nweight", option = "inferno", breaks = c(0, 0.5, 1)) +
    scale_x_continuous("reduced dimension 1") +
    scale_y_continuous("reduced dimension 2") +
    ggtitle(paste("Pattern", i)) +
    mytheme
  plot_list[[i]] = p
}

pdf(paste0("HSCR_proj_Linnarsson_multipage.pdf"), width = 7.25, height = 2.5, useDingbats = F)
for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
}
dev.off()

pdf("HSCR_proj_Linnarsson_singlepage.pdf", width = 30, height = 30, useDingbats = F)
do.call("grid.arrange", c(plot_list, ncol = 5))
dev.off()


pdf("HSCR_proj_Linnarsson_pg1.pdf", width = 12, height = 18, useDingbats = F)
do.call("grid.arrange", c(plot_list[1:18], ncol = 2))
dev.off()

pdf("HSCR_proj_Linnarsson_pg2.pdf", width = 12, height = 18, useDingbats = F)
do.call("grid.arrange", c(plot_list[19:36], ncol = 2))
dev.off()

pdf("HSCR_proj_Linnarsson_pg3.pdf", width = 12, height = 18, useDingbats = F)
do.call("grid.arrange", c(plot_list[37:54], ncol = 2))
dev.off()

pdf("HSCR_proj_Linnarsson_pg4.pdf", width = 12, height = 16, useDingbats = F)
do.call("grid.arrange", c(plot_list[55:70], ncol = 2))
dev.off()
