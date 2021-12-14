#########################################################################################################################
# HSCR GWAS genes
#########################################################################################################################


HSCR_GWAS_genes <- c("Ednrb", "Gfra1", "Ikbkap", "Sox10", "Zeb2", 
                "Fam213a", "Kif1bp", "Phox2b", "Tcf4", "Ubr4", 
                "Ece1", "L1cam", "Nrg1", "Sema3c", "Sema3d", 
                "Acss2", "Ret", "Adamts17", "Eno3", "Sh3pxd2a", 
                "Slc27a4", "Edn3", "Gdnf", "Nrtn")

HSCR_GWAS_genes <- sort(HSCR_GWAS_genes)

HSCR_GWAS_genes <- HSCR_GWAS_genes[HSCR_GWAS_genes %in% lookupGeneName(dat.filtered, expressed_genes)]

violindat <- merge(pData(dat.filtered)[,c("genotype", "cluster")],
                    t(exprs(dat.filtered)[rownames(fData(dat.filtered)[fData(dat.filtered)$gene_short_name %in% HSCR_GWAS_genes,]),]),
                   by = 0)

rownames(violindat) <- violindat$Row.names
violindat <- violindat[,-1]

melted_violindat <- melt(violindat, id.vars = c("genotype", "cluster"))

plot_list = list()

for(i in seq_along(HSCR_GWAS_genes)){
  p = melted_violindat[melted_violindat$variable == lookupGeneId(dat.filtered, HSCR_GWAS_genes)[i],] %>%
    ggplot(., aes(x = cluster, y = log10(value + 1), fill = cluster)) +
    geom_violin(trim = T, scale = "width") +
    stat_summary(fun.y = mean, geom = "point", size = 1) +
    scale_fill_manual(values = c(brewer.pal(8, "Set2"), brewer.pal(5, "Set3")[5])) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = "none", 
          plot.title = element_text(face = "plain", size = 12),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    ) +
    ggtitle(HSCR_GWAS_genes[i])
  plot_list[[i]] = p
}

pdf("SF6C-HSCR_GWAS_genes_expressed_violin_log=T_trim=T.pdf", height = 5, width = 14, useDingbats = F)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
             plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
             plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]],
             plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]],
             plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]],
             plot_list[[21]], plot_list[[22]], ncol = 8
)
dev.off()



plot_list = list()


for(i in seq_along(HSCR_GWAS_genes)){
  p = melted_violindat[melted_violindat$variable == lookupGeneId(dat.filtered, HSCR_GWAS_genes)[i],] %>%
    ggplot(., aes(x = genotype, y = log10(value + 1), fill = genotype)) +
    geom_violin(trim = T, scale = "width") +
    stat_summary(fun.y = mean, geom = "point", size = 1) +
    scale_fill_manual(values = c("gray70", "darkred")) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = "none",
          aspect.ratio = 1,
          plot.title = element_text(face = "plain", size = 12),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    ) +
    ggtitle(HSCR_GWAS_genes[i])
  plot_list[[i]] = p
}

pdf("SF6B-HSCR_GWAS_genes_expressed_violin_log=T_trim=T_genotype.pdf", height = 7, width = 12, useDingbats = F)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
             plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
             plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]],
             plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]],
             plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]],
             plot_list[[21]], plot_list[[22]], ncol = 6
)
dev.off()


pdf("SF6A-RetGRN-UMAPs.pdf", height = 7, width = 9, useDingbats = F)
plotUMAP(dat.filtered, markers = HSCR_GWAS_genes, scaled = T, size = 0.2, sort = F, ncol = 6) +
  scale_color_viridis(option = "inferno", breaks = c(0,1)) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, vjust = -2),
        legend.direction = "horizontal",
        legend.position = c(0.7, 0.1),
        legend.key.height = unit(1, "line"),
        legend.key.width = unit(1.5, "line"))
dev.off()

#########################################################################################################################
# Ret GRN genes
#########################################################################################################################

GRN_genes <- c("Ret", "Ednrb", "Sox10", "Gata2", "Rarb", "Nkx2-5", "Gfra1", "Gdnf", "Cbl")

GRN_genes <- sort(GRN_genes)

GRN_genes <- GRN_genes[GRN_genes %in% lookupGeneName(dat.filtered, expressed_genes)]
GRN_genes <- GRN_genes[!GRN_genes %in% HSCR_GWAS_genes]

violindat2 <- merge(pData(dat.filtered)[,c("genotype", "cluster")],
                   t(exprs(dat.filtered)[rownames(fData(dat.filtered)[fData(dat.filtered)$gene_short_name %in% GRN_genes,]),]),
                   by = 0)

rownames(violindat2) <- violindat2$Row.names
violindat2 <- violindat2[,-1]

melted_violindat2 <- melt(violindat2, id.vars = c("genotype", "cluster"))

plot_list = list()

for(i in seq_along(GRN_genes)){
  p = melted_violindat2[melted_violindat2$variable == lookupGeneId(dat.filtered, GRN_genes)[i],] %>%
    ggplot(., aes(x = cluster, y = log10(value + 1), fill = cluster)) +
    geom_violin(trim = T, scale = "width") +
    stat_summary(fun.y = mean, geom = "point", size = 1) +
    scale_fill_manual(values = c(brewer.pal(8, "Set2"), brewer.pal(5, "Set3")[5])) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = "none", 
          plot.title = element_text(face = "plain", size = 12),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    ) +
    ggtitle(GRN_genes[i])
  plot_list[[i]] = p
}

pdf("SF6F-GRN_genes_expressed_violin_log=T_trim=T.pdf", height = 1.7, width = 5.3, useDingbats = F)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], nrow = 1
)
dev.off()



plot_list = list()


for(i in seq_along(GRN_genes)){
  p = melted_violindat2[melted_violindat2$variable == lookupGeneId(dat.filtered, GRN_genes)[i],] %>%
    ggplot(., aes(x = genotype, y = log10(value + 1), fill = genotype)) +
    geom_violin(trim = T, scale = "width") +
    stat_summary(fun.y = mean, geom = "point", size = 1) +
    scale_fill_manual(values = c("gray70", "darkred")) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = "none",
          aspect.ratio = 1,
          plot.title = element_text(face = "plain", size = 12),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    ) +
    ggtitle(GRN_genes[i])
  plot_list[[i]] = p
}

pdf("SF6E-GRN_genes_expressed_violin_log=T_trim=T_genotype.pdf", height = 2, width = 5.3, useDingbats = F)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], nrow = 1
)
dev.off()


pdf("SF6D-RetGRN-UMAPs.pdf", height = 2, width = 6, useDingbats = F)
plotUMAP(dat.filtered, markers = GRN_genes, scaled = T, size = 0.2, sort = F, ncol = 6) +
  scale_color_viridis(option = "inferno", breaks = c(0,1)) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
dev.off()
