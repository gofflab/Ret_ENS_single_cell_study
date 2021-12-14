pdf("Read_distribution_by_condition.pdf", height = 3, width = 8, useDingbats = F)
pData(dat.relative) %>%
  mutate(group = paste(age, sex)) %>%
  ggplot(.) + 
  geom_density(aes(x = log10(total_reads), fill = genotype), alpha = 0.4) +
  facet_wrap(~group, nrow = 1) +
  scale_fill_manual(breaks = c("het", "hom"), values = c("gray10", "darkred")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  ggtitle("Total Reads Distributions")
dev.off()

##################################################################################

pdf("Read_distribution_by_condition_filtered.pdf", height = 3, width = 8, useDingbats = F)
pData(dat.filtered) %>%
  mutate(group = paste(age, sex)) %>%
  ggplot(.) + 
  geom_density(aes(x = log10(total_reads), fill = genotype), alpha = 0.4) +
  facet_wrap(~group, nrow = 1) +
  scale_fill_manual(breaks = c("het", "hom"), values = c("gray10", "darkred")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  ggtitle("Total Reads Distributions")
dev.off()

##################################################################################

pdf("PC1vsPC2-batch.pdf", height = 4, width = 5, useDingbats = F)
ggbiplot(dat.filtered.pca, var.axes = F, groups = pData(dat.filtered)$batch, alpha = 0, ellipse = T) +
  ggtitle("PC1 vs. PC2") + 
  scale_color_viridis("Batch", discrete = T, option = "B") +  
  geom_point(aes(color = pData(dat.filtered)$batch), size = 0.6) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.7, 0.2),
        plot.title = element_text(size = 12), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))
dev.off()

#################################################################################


mytable <- melt(table(pData(dat.relative)[!rownames(pData(dat.relative)) %in% rownames(pData(dat.filtered)),c("age", "sex", "genotype")]))
mytable$all <- melt(table(pData(dat.relative)[,c("age", "sex", "genotype")]))$value
mytable$removed <- 1

mytable2 <- melt(table(pData(dat.filtered)[,c("age", "sex", "genotype")]))
mytable2$all <- melt(table(pData(dat.relative)[,c("age", "sex", "genotype")]))$value
mytable2$removed <- 0

mytable3 <- rbind(mytable2, mytable)
mytable3$age <- as.factor(mytable3$age)
mytable3$sex <- as.factor(mytable3$sex)
mytable3$genotype <- as.factor(mytable3$genotype)

pdf("Cell_distribution_filtered.pdf", height = 4, width = 8, useDingbats = F)
ggplot(mytable3, aes(x = genotype, y = value)) + 
  geom_col(aes(fill = genotype, alpha = -removed)) + 
  facet_wrap(~age + sex, nrow = 1) +
  scale_y_continuous("# Cells", expand = expand_scale(mult = c(0, .05))) +
  scale_fill_manual(values = c("gray70", "darkred")) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  geom_text(stat = "identity", aes(y = all, label = all), vjust = -0.4) +
  geom_text(stat = "identity", aes(y = value, label = value), vjust = -0.4) +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(b = 5)))
dev.off()

#############################################################


my_plot_ordering_genes <- function(disp_table, genes){
  p <- ggplot(disp_table, aes(x = log10(mean_expression), 
                              y = log10(dispersion_empirical))) + 
  geom_point(aes(color = ifelse(disp_table$gene_id %in% genes, "1", "0")),
             size = 0.2, alpha = 0.2) +
  geom_line(aes(y = log10(dispersion_fit)), color = "red") +
  scale_color_manual("", values = c("gray10", "red")) +
  theme(aspect.ratio = 1,
        legend.position = "none")
  p
}

###################################################################

ggbiplot(pca.2, var.axes = F, groups = pData(dat.filtered)$batch, alpha = 0, ellipse = T) +
  ggtitle("PC1 vs. PC2") + 
  scale_color_viridis("Batch", discrete = T, option = "B") +  
  geom_point(aes(color = pData(dat.filtered)$batch), size = 0.6) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.7, 0.2),
        plot.title = element_text(size = 12), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))


ggbiplot(pca.3, var.axes = F, groups = pData(dat.filtered)$batch, alpha = 0, ellipse = T) +
  ggtitle("PC1 vs. PC2") + 
  scale_color_viridis("Batch", discrete = T, option = "B") +  
  geom_point(aes(color = pData(dat.filtered)$batch), size = 0.6) +
#  scale_x_continuous(limits = c(-0.5, 0.5)) +
#  scale_y_continuous(limits = c(-0.5, 0.5)) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.7, 0.2),
        plot.title = element_text(size = 12), 
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))

########################################################################

pdf("UMAP_batch.pdf", width = 3, height = 2.5, useDingbats = F)
plotUMAP(dat.filtered, color = "batch", size = 0.1) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 3))) + 
  scale_color_viridis("Batch", discrete = T, option = "B") +
  theme(aspect.ratio = 1, 
        legend.position = c(0.5, 0.2),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.key.height = unit(0.8, "line"),
        legend.key.width = unit(0.5, "line")
  )
dev.off()

pdf("UMAP_age.pdf", width = 3, height = 2.5, useDingbats = F)
plotUMAP(dat.filtered, color = "age", size = 0.1) +
  guides(color = guide_legend(title = "Age", override.aes = list(size = 3))) + 
  scale_color_manual(values = brewer.pal(4, "Set1")[3:4]) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.7, 0.2),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.key.height = unit(0.8, "line"),
        legend.key.width = unit(0.5, "line")
  )
dev.off()
