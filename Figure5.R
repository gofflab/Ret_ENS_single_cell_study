glia_terms <- c("GO:0042063", "GO:0045165", "GO:0014033", "GO:0061351", "GO:0045787", 
                "GO:0042552", "GO:0007088", "GO:0045786", "GO:0000132", "GO:0051302", 
                "GO:0000082")

ego_glia_subset <- ego_glia_genotype@result[glia_terms,]

ego_glia_subset$Description <- factor(ego_glia_subset$Description, levels = rev(ego_glia_subset$Description))

pdf("F4B-DEG_glia!cluster9_genotype_GO_plot.pdf", height = 2.5, width = 6, useDingbats = F)
ggplot(ego_glia_subset) +
  geom_col(aes(x = Description, y = Count, fill = qvalue)) +
  scale_x_discrete("GO:BP Term") +
  scale_y_continuous("# Genes", expand = c(0,0)) +
  scale_fill_viridis("q-value", limits = c(0, 0.05), breaks = c(0, 0.05), option = "plasma") +
  coord_flip() +
  theme(legend.position = c(0.77, 0.35), 
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
        axis.title.y = element_text(margin = margin(r = -10), hjust = 0.6),
        axis.title.x = element_text(margin = margin(t = -2)))
dev.off()

heatmap_data <- log10(exprs(dat.filtered)[glia_genotype_diff_test_sig_genes, pData(dat.filtered)$celltype == "progenitor/glia" & !pData(dat.filtered)$cluster == 9] + 1) # log transform
heatmap_data <- t(apply(heatmap_data, 1, function(x){ # Normalize gene expr by dividing by max
  x/max(x)
}))

heatmap_data <- heatmap_data[order.dendrogram(as.dendrogram(hclust(dist(heatmap_data)))),order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data)))))]
rownames(heatmap_data) <- lookupGeneName(dat.filtered, rownames(heatmap_data))

gene_ids <- str_split(ego_glia_genotype@result[glia_terms, "geneID"], "\\/")
names(gene_ids) <- ego_glia_genotype@result[glia_terms, "Description"]

tmp <- lapply(gene_ids, function(x){lookupGeneName(dat.filtered, glia_genotype_diff_test_sig_genes) %in% x})
heatmap_annotation_row <- do.call(cbind, tmp)
heatmap_annotation_row <- as.data.frame(ifelse(heatmap_annotation_row, 1, 0))

heatmap_annotation_row <- heatmap_annotation_row[,order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_annotation_row)))))]

CC_regulation_terms <- c("GO:0045787", #positive regulation of cell cycle
             "GO:0045786", #negative regulation of cell cycle
             "GO:0007088", #regulation of mitotic nuclear division
             "GO:0000082", #G1/S transition of mitotic cell cycle
             "GO:0000132", #establishment of mitotic spindle orientation
             "GO:0051302" #regulation of cell division
)

heatmap_annotation_row$CC_regulation <- apply(heatmap_annotation_row[,ego_glia_subset[CC_regulation_terms, "Description"]], 1, function(x){
  ifelse(sum(x) > 0, 1, 0)
})

heatmap_annotation_row <- as.data.frame(apply(heatmap_annotation_row, 2, as.factor))
rownames(heatmap_annotation_row) <- lookupGeneName(dat.filtered, glia_genotype_diff_test_sig_genes)

heatmap_annotation_col <- pData(dat.filtered)[,c("genotype", "age", "cluster")]

writeLines(paste0("\"", 1:9, "\" = \"", c(brewer.pal(8, "Set2"), brewer.pal(5, "Set3")[5]), "\""), sep = ",\n")

heatmap_colors <- list(
  "genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"),
  "cluster" = c("1" = "#66C2A5",
                "2" = "#FC8D62",
                "3" = "#8DA0CB",
                "4" = "#E78AC3",
                "5" = "#A6D854",
                "6" = "#FFD92F",
                "7" = "#E5C494",
                "8" = "#B3B3B3",
                "9" = "#80B1D3"),
  "age" = c("E12.5" = "#4daf4a",
            "E14.5" = "#984ea3"),
  
  "CC_regulation" = c("1" = "#E7298A", "0" = "grey90")
  )

#colnames(heatmap_annotation_col) <- c("Cell Type","Genotype")

pdf("DEG_glia!cluster9_genotype_heatmap_test.pdf", height = 20, width = 20)
pheatmap(mat = heatmap_data,
         scale="none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         drop_levels = T,
         cluster_rows = F,
         cluster_cols = F,
         color = viridis(100),
         annotation_col = heatmap_annotation_col,
         annotation_row = heatmap_annotation_row[,"CC_regulation", drop = F],
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = T
)
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

  
