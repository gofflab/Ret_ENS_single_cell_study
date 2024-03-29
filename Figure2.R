#Figure 2A
heatmap_data <- log10(exprs(dat.filtered)[genotype_diff_test_sig_genes,] + 1) # log transform
heatmap_data <- t(apply(heatmap_data, 1, function(x){ # Normalize gene expr by dividing by max
  x/max(x)
}))

rownames(heatmap_data) <- lookupGeneName(dat.filtered, rownames(heatmap_data))

heatmap_data <- heatmap_data[order.dendrogram(as.dendrogram(hclust(dist(heatmap_data)))),order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data)))))]

heatmap_annotation <- pData(dat.filtered)[, "genotype", drop = F]
colnames(heatmap_annotation) <- "Genotype"

heatmap_colors <- list(
  "Genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"))

pdf("DEG_all_cells_genotype_heatmap.pdf", height = 24, width = 24)
pheatmap(mat = heatmap_data,
         scale="none",
         show_rownames = TRUE,
         labels_row = unlist(lapply(rownames(heatmap_data), function(x){ifelse(x %in% unlist(gene_ids), x, "")})),
         show_colnames = FALSE,
         drop_levels = FALSE,
         cluster_rows = T,
         cluster_cols = T,
         color = viridis(100),
         annotation_col = heatmap_annotation,
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE,
         fontsize_row = 2,
         treeheight_row = 20,
         treeheight_col = 20
)
dev.off()

heatmap_data_het <- heatmap_data[,rownames(pData(dat.filtered)[rownames(pData(dat.filtered)) %in% colnames(heatmap_data) & pData(dat.filtered)$genotype == "het",])]
heatmap_data_het <- heatmap_data_het[,order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data_het)))))]

heatmap_data_hom <- heatmap_data[,rownames(pData(dat.filtered)[rownames(pData(dat.filtered)) %in% colnames(heatmap_data) & pData(dat.filtered)$genotype == "hom",])]
heatmap_data_hom <- heatmap_data_hom[,order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data_hom)))))]

heatmap_data <- cbind(heatmap_data_hom, heatmap_data_het)

heatmap_data <- heatmap_data[order.dendrogram(as.dendrogram(hclust(dist(heatmap_data)))),]


gene_ids <- str_split(ego_genotype_subset$geneID, "\\/")
names(gene_ids) <- ego_genotype_subset$Description

tmp <- lapply(gene_ids, function(x){lookupGeneName(dat.filtered, genotype_diff_test_sig_genes) %in% x})
heatmap_annotation_row <- do.call(cbind, tmp)
heatmap_annotation_row <- as.data.frame(ifelse(heatmap_annotation_row, 1, 0))

neurotransmitters <- c("GO:0001505", #regulation of neurotransmitter levels
                            "GO:0006584" #catecholamine metabolic process
)

axonogenesis <- c("GO:0061564", #axon development
                  "GO:0050808", #synapse organization
                  "GO:1990138" #neuron projection extension
                  )

cellcycle <- c("GO:0045787", #positive regulation of cell cycle
               "GO:0051302", #regulation of cell division
               "GO:0007088" #regulation of mitotic nuclear division
               )

heatmap_annotation_row$neurotransmitters <- apply(heatmap_annotation_row[,ego_genotype_subset[neurotransmitters, "Description"]], 1, function(x){
  ifelse(sum(x) > 0, 1, 0)
})

heatmap_annotation_row$axonogenesis <- apply(heatmap_annotation_row[,ego_genotype_subset[axonogenesis, "Description"]], 1, function(x){
  ifelse(sum(x) > 0, 1, 0)
})

heatmap_annotation_row$cellcycle <- apply(heatmap_annotation_row[,ego_genotype_subset[cellcycle, "Description"]], 1, function(x){
  ifelse(sum(x) > 0, 1, 0)
})


heatmap_annotation_row <- heatmap_annotation_row[,c("cellcycle", "positive regulation of neuron differentiation", "anterior/posterior pattern specification")]


heatmap_annotation_row <- heatmap_annotation_row[,order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_annotation_row)))))]

heatmap_annotation_row <- as.data.frame(apply(heatmap_annotation_row, 2, as.factor))
rownames(heatmap_annotation_row) <- lookupGeneName(dat.filtered, genotype_diff_test_sig_genes)

heatmap_annotation_col <- pData(dat.filtered)[,c("genotype", "celltype")]

writeLines(paste0("\"", colnames(heatmap_annotation_row), "\" = c(\"1\" = \"black\", \"0\" = \"grey90\")"), sep = ",\n")

heatmap_colors <- list(
  "genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"),
  "celltype"  =  c("neuron" = "#E69F00", "progenitor/glia" = "#56B4E9"),
    
  "cellcycle" = c("1" = "#E7298A", "0" = "grey90"),
  "positive regulation of neuron differentiation" = c("1" = "#66a61e", "0" = "grey90"),
  "anterior/posterior pattern specification" = c("1" = "#a6761d", "0" = "grey90")
)

pdf("F2A-DEG_all_cells_genotype_heatmap_annotated.pdf", height = 12, width = 12)
pheatmap(mat = heatmap_data,
         scale="none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         drop_levels = T,
         cluster_rows = F,
         cluster_cols = T,
         color = viridis(100),
         annotation_col = heatmap_annotation_col,
         annotation_row = heatmap_annotation_row,
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = F,
         treeheight_col = 30
)
dev.off()


#Figure 2B

egENSEMBL <- toTable(org.Mm.egENSEMBL)

genotype_diff_test_entrez <- gsub("\\..*","", genotype_diff_test_sig_genes)
m <- match(genotype_diff_test_entrez, egENSEMBL$ensembl_id)
genotype_diff_test_entrez <- egENSEMBL$gene_id[m]
ego_genotype <- enrichGO(gene = genotype_diff_test_entrez,
                         OrgDb =       'org.Mm.eg.db',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 0.05, 
                         readable      = TRUE)

ego_genotype_subset <- ego_genotype@result[ego_genotype@result$ID %in% c("GO:0007219", "GO:0042063", "GO:0045165", "GO:0061564", "GO:0050768", "GO:0045666", "GO:0001505", "GO:0006584", "GO:0061351", "GO:0014033", "GO:0045787", "GO:0030335", "GO:0070371", "GO:0042552", "GO:0009952", "GO:0010718", "GO:0098727", "GO:0043405", "GO:0051302", "GO:0050808", "GO:1990138", "GO:0007088"),]

ego_genotype_subset$Description <- factor(ego_genotype_subset$Description, levels = rev(ego_genotype_subset$Description))

pdf("DEG_all_cells_genotype_GO_plot.pdf", height = 4, width = 8, useDingbats = F)
ggplot(ego_genotype_subset) +
  geom_col(aes(x = Description, y = Count, fill = qvalue)) +
  scale_x_discrete("GO:BP Term") +
  scale_y_continuous("# Genes", expand = c(0,0)) +
  scale_fill_viridis("q-value", limits = c(0, 0.05), breaks = c(0, 0.05), option = "plasma") +
  coord_flip() +
  theme(legend.position = c(0.7, 0.2), 
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
        axis.title.y = element_text(margin = margin(r = -10)),
        axis.title.x = element_text(margin = margin(t = -2)))
dev.off()


#Figure 2C

antposDEG <- unlist(str_split(ego_genotype_subset[ego_genotype_subset$ID == "GO:0009952", "geneID"], "\\/"))
#GO:0009952 = anterior/posterior pattern specification
names(antposDEG) <- lookupGeneId(dat.filtered, antposDEG)

hoxDEG <- c("Hoxa2", "Hoxb2", "Hoxd3", "Hoxa4", "Hoxa5", "Hoxb5", "Hoxb6")
names(hoxDEG) <- lookupGeneId(dat.filtered, hoxDEG)

violindat <- merge(pData(dat.filtered)[,c("genotype", "cluster")],
#                   t(exprs(dat.filtered)[names(hoxDEG),]),
                  t(exprs(dat.filtered)[names(antposDEG),]),
                  by = 0)

rownames(violindat) <- violindat$Row.names
violindat <- violindat[,-1]

melted_violindat <- melt(violindat, id.vars = c("genotype", "cluster"))

plot_list = list()

#for(i in seq_along(antposDEG)){
  p = melted_violindat[melted_violindat$variable == names(antposDEG)[i],] %>%
    ggplot(., aes(x = genotype, y = log10(value + 1), fill = genotype)) +
    geom_violin(trim = T, scale = "width") +
    #geom_jitter() +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none", 
          legend.title = element_text(face = "plain"),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(antposDEG[i])
  plot_list[[i]] = p
}


pdf("F2C-antposDEG-by-genotype.pdf", width = 8, height = 12)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
             plot_list[[4]], plot_list[[5]], plot_list[[6]], 
             plot_list[[7]], plot_list[[8]], plot_list[[9]], 
             plot_list[[10]], plot_list[[11]], plot_list[[12]], 
             plot_list[[13]], plot_list[[14]], ncol = 3)
dev.off()

