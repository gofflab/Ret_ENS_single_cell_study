plotUMAP(dat.filtered, color = "CellType") +
  scale_color_manual("Cell Type", values = c(cbPalette[c(1:3, 6)], "#999999")) +
  theme(aspect.ratio = 1, legend.position = c(0.8, 0.2))
plotUMAP(tmp, color = "subset_pred_description") + 
  scale_color_manual("Cell Type", values = cbPalette[c(1, 5, 6, 7)]) +
  theme(aspect.ratio = 1, , legend.position = c(0.7, 0.2))

plotUMAP(dat.filtered, color = "cluster") +
#  scale_color_manual("Cluster", values = colorRampPalette(brewer.pal(9, "Set1"))(10), guide = guide_legend(ncol = 2)) +
  scale_color_brewer("Cluster", palette = "Set1", guide = guide_legend(ncol = 2)) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = c(0.8, 0.2))

tmp <- cluster_specificity[cluster_specificity$cluster == 4, ]
ggplot(tmp, aes(x = log10(mean_expr + 1), y = specificity, label = ifelse(tmp$gene_short_name %in% c("Tac1", "Ecel1", "Scng", "Th", "Chgb", "Mapt", "Prph", "Stmn3", "Sncg", "Snap25", "Meg3", "Stmn2", "Rtn1"), tmp$gene_short_name, ""))) +
  geom_point() +
  geom_smooth() +
  geom_text(nudge_y = 0.015) +
  ggtitle(paste("Specificity for Cluster 4")) +
  theme(aspect.ratio = 1)

tmp <- cluster_specificity[cluster_specificity$cluster == 9, ]
ggplot(tmp, aes(x = log10(mean_expr + 1), y = specificity, label = ifelse(tmp$gene_short_name %in% c("Nos1", "Ntng1", "Sncg", "Mapt", "Prph", "Vip", "Cartpt", "Myl1", "Etv1"), tmp$gene_short_name, ""))) +
  geom_point() +
  geom_smooth() +
  geom_text(nudge_y = 0.015) +
  ggtitle(paste("Specificity for Cluster 9")) +
  theme(aspect.ratio = 1)

tmp <- cluster_specificity[cluster_specificity$cluster == 10, ]
ggplot(tmp, aes(x = log10(mean_expr + 1), y = specificity, label = ifelse(tmp$gene_short_name %in% c("Plp1", "Myl9", "Fabp7", "Postn", "Sox5", "Gja1", "Tagln", "Pls3", "Gpc3", "Sparc", "Ccnd2", "Matn4"), tmp$gene_short_name, ""))) +
  geom_point() +
  geom_smooth() +
  geom_text(nudge_y = 0.015) +
  ggtitle(paste("Specificity for Cluster 10")) +
  theme(aspect.ratio = 1)

ggplot(pData(dat.filtered)) +
  geom_bar(aes(x = cluster, fill = genotype), position = "dodge") +
  scale_fill_brewer("Genotype", palette = "Set1") +
  labs(x = "Cluster", y = "Count")

ggplot(pData(dat.filtered)[pData(dat.filtered)$cluster %in% c(4, 9),]) +
  geom_bar(aes(x = cluster, fill = genotype), position = "dodge") +
  scale_fill_brewer("Genotype", palette = "Set1") +
  labs(x = "Cluster", y = "Count")

plotUMAP(dat.filtered, color = "genotype") +
  theme(aspect.ratio = 1, legend.position = c(0.8, 0.2))



tmp <- t(scale(t(exprs(dat.filtered[genotypeCellType_diff_test_sorted_sig_genes,])),scale = T,center = T))
heatmap_data <- log(tmp -min(tmp) + 1)

myLower <- 0
myUpper <- 2.5

heatmap_data[heatmap_data < myLower] <- myLower
heatmap_data[heatmap_data > myUpper] <- myUpper

rownames(heatmap_data) <- lookupGeneName(dat.filtered,rownames(heatmap_data))

heatmap_annotation <- pData(dat.filtered)[,c("CellType","genotype")]

heatmap_colors <- list(
  "Cell Type" =  c("Acetylcholinergic Neuron" = "#E69F00", "Neuron" = "#56B4E9", "Noradrenergic Neuron" = "#009E73", "Progenitor/Glia" = "#D55E00", "Transition Cell" = "#999999"),
  "Genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"))

colnames(heatmap_annotation) <- c("Cell Type","Genotype")

pheatmap(mat = heatmap_data,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = FALSE,
         drop_levels = FALSE,
         breaks = seq(myLower, myUpper, length = 100),
         annotation_col = heatmap_annotation,
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         fontsize_row = 8,
         cutree_rows = 7,
         cutree_cols = 7
)



#tmp <- t(scale(t(exprs(dat.filtered[genotype_diff_test_sorted_sig_genes,])), scale = T, center = T))

tmp <- t(scale(t(exprs(dat.filtered[rownames(genotype_diff_test_sorted[genotype_diff_test_sorted$qval < 0.001,]),])), scale = T, center = T))
heatmap_data <- log10(tmp - min(tmp) + 1)

rownames(heatmap_data) <- lookupGeneName(dat.filtered, rownames(heatmap_data))

#heatmap_data <- heatmap_data[order.dendrogram(as.dendrogram(hclust(dist(heatmap_data)))),order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data)))))]

myLower <- 0.1
myUpper <- 0.9

heatmap_data[heatmap_data < myLower] <- myLower
heatmap_data[heatmap_data > myUpper] <- myUpper

heatmap_annotation <- pData(dat.filtered)[, c("CellType", "genotype")]

heatmap_colors <- list(
  "Cell Type" =  c("Acetylcholinergic Neuron" = "#E69F00", "Neuron" = "#56B4E9", "Noradrenergic Neuron" = "#009E73", "Progenitor/Glia" = "#D55E00", "Transition Cell" = "#999999"),
  "Genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"))

colnames(heatmap_annotation)<-c("Cell Type","Genotype")

#pdf("Heatmap_all_cells_wrt_genotype.pdf",height=25,width=18)
#png("Heatmap_all_cells_wrt_genotype.png",height=15,width=20, units = "in", res = 400)
pheatmap(mat = heatmap_data,
         scale="none",
         show_rownames = TRUE,
         show_colnames = FALSE,
         drop_levels = FALSE,
         #cluster_rows = F,
         #cluster_cols = F,
         breaks=seq(myLower, myUpper, length = 100),
         annotation_col = heatmap_annotation,
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         fontsize_row = 8,
         treeheight_col = 0,
         treeheight_row = 0,
         #cutree_rows = 11,
         cutree_cols = 10
)





heatmap_annotation <- pData(dat.glia)[, c("genotype", "cluster")]
heatmap_colors <- list(
  "Genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"),
  "Cluster" = c("1" = "#E41A1C", "2" = "#596A98", "3" = "#449B75", "4" = "#6B886D", "5" = "#AC5782", "6" = "#FF7F00", "7" = "#FFE528", "8" = "#C9992C", "9" = "#C66764", "10" = "#E485B7", "11" = "#999999")
)

colnames(heatmap_annotation) <- c("Genotype","Cluster")
heatmap_annotation <- heatmap_annotation %>%
  select("Genotype","Cluster") #order the columns

tmp <- t(scale(t(exprs(dat.glia[glia_genotype_diff_test_sigGenes,])), scale = T, center = T))
heatmap_data <- log(tmp -min(tmp) + 1)
#heatmap_data<-tmp -min(tmp)
#myLower<-(0)
#myUpper<-(3)
#heatmap_data[heatmap_data<myLower]<-myLower
#heatmap_data[heatmap_data>myUpper]<-myUpper
rownames(heatmap_data) <- lookupGeneName(dat.glia, rownames(heatmap_data))

#pdf("posterFig_Heatmap_glia_wrt_genotype.pdf",height=36,width=17)
pheatmap(mat = heatmap_data,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = FALSE,
         drop_levels = FALSE,
         #breaks = seq(myLower, myUpper, length = 100),
         annotation_col = heatmap_annotation,
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         border_color = NA,
         fontsize = 4,
         cutree_rows = 5,
         cutree_cols = 7,
         treeheight_row = 0,
         treeheight_col = 0
)

plotUMAP(dat.filtered, markers = "Ret") +
  theme(aspect.ratio = 1, legend.position = c(0.8, 0.2))

bardat <- as.data.frame(table(pData(dat.filtered)[,c("cluster", "genotype")]))
bardat$proportion <- NA
for(i in 1:8){
  bardat[bardat$cluster == i, "proportion"] <- bardat[bardat$cluster == i, "Freq"]/sum(bardat[bardat$cluster == i, "Freq"])
}

ggplot(bardat) +
  geom_col(aes(x = cluster, y = proportion, fill = genotype)) +
  geom_hline(yintercept = table(pData(dat.filtered)$genotype)[2]/nrow(pData(dat.filtered)), lty = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Cluster", y = "Proportion", fill = "Genotype") +
  scale_fill_brewer(palette = "Set1")

violindat1 <- merge(pData(dat.filtered)[,c("genotype", "cluster")], t(exprs(dat.filtered)[lookupGeneId(dat.filtered, "Ret"),,drop = F]), by = 0)
rownames(violindat1) <- violindat1$Row.names
violindat1 <- violindat1[,-1]
names(violindat1) <- c("genotype", "cluster", "expression")
violindat1$gene_name <- "Ret"

violindat2 <-merge(pData(dat.filtered)[,c("genotype", "cluster")], t(exprs(dat.filtered)[lookupGeneId(dat.filtered, "tCFP"),,drop = F]), by = 0)
rownames(violindat2) <- violindat2$Row.names
violindat2 <- violindat2[,-1]
names(violindat2) <- c("genotype", "cluster", "expression")
violindat2$gene_name <- "CFP"

violindat <- rbind(violindat1, violindat2)

ggplot(violindat1, aes(x = cluster, y = expression, color = genotype)) +
#  geom_boxplot(aes(x = cluster, y = log10(violindat$expression + 1), color = genotype)) +
  geom_point(position = position_jitterdodge()) +
  geom_boxplot(position = position_dodge()) +
#  facet_wrap(~gene_name) +
  scale_y_continuous(expand = c(0, 0))

ggplot(violindat, aes(x = cluster, y = expression, fill = genotype)) +
#  geom_col(aes(x = cluster, y = log10(violindat$expression + 1), fill = genotype), position = "dodge") +
  geom_col(position = "dodge") +
  facet_wrap(~gene_name) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1")

ggplot(violindat[violindat$gene_name == "Ret",]) +
  geom_violin(aes(x = cluster, y = log10(violindat[violindat$gene_name == "Ret","expression"] + 1), color = genotype, fill = genotype)) +
#  geom_violin(aes(x = cluster, y = expression, color = genotype, fill = genotype)) +
#  facet_wrap(~gene_name) +
  scale_y_continuous(expand = c(0, 0))
                     
plotUMAP(dat.filtered, markers = "Slc27a4", size = 1) + 
  scale_color_gradient(low = "gray", high = "darkorange1") + 
  theme_classic() + 
  theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())