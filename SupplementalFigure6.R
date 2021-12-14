ego_genotype_celltype_subset1 <- ego_genotype_celltype@result[ego_genotype_celltype@result$ID %in% 
                                                                c("GO:0007059", "GO:0006260", "GO:0045787", "GO:0044770", "GO:0006281", 
                                                                  "GO:0045786", "GO:1902850", "GO:0007051", "GO:0051383", "GO:0071103", 
                                                                  "GO:0051321", "GO:0033044", "GO:0051656", "GO:0006310", "GO:0034502", 
                                                                  "GO:0000910", "GO:0051302", "GO:0007219", "GO:0042136", "GO:0010718",
                                                                  "GO:0042063"),]

ego_genotype_celltype_subset1$Description <- factor(ego_genotype_celltype_subset1$Description, levels = rev(ego_genotype_celltype_subset1$Description))

pdf("SF6B-DEG_all_cells_genotypeXcelltype_GO_plot.pdf", height = 3.5, width = 8, useDingbats = F)
ggplot(ego_genotype_celltype_subset1) +
  geom_col(aes(x = Description, y = Count, fill = qvalue)) +
  scale_x_discrete("GO:BP Term") +
  scale_y_continuous("# Genes", expand = c(0,0)) +
  scale_fill_viridis("q-value", limits = c(0, 0.05), breaks = c(0, 0.05), option = "plasma") +
  coord_flip() +
  theme(legend.position = c(0.7, 0.25), 
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
        #        axis.title.y = element_text(margin = margin(r = -10)),
        axis.title.x = element_text(margin = margin(t = -2)))
dev.off()

write.csv(ego_genotype_celltype@result[ego_genotype_celltype@result$qvalue < 0.05,], "DEG_all_cells_genotypeXcelltype_GO.csv")

gene_ids <- str_split(ego_genotype_celltype_subset1$geneID, "\\/")
names(gene_ids) <- str_remove(ego_genotype_celltype_subset1$ID, ":")
tmp <- lapply(gene_ids, function(x){lookupGeneName(dat.filtered, genotypeCellType_diff_test_sig_genes) %in% x})
heatmap_annotation_row <- do.call(cbind, tmp)
heatmap_annotation_row <- as.data.frame(ifelse(heatmap_annotation_row, 1, 0))

S_terms <- c("GO0006260", #DNA replication
             "GO0006281" #DNA repair
)
EarlyMidM_terms <- c("GO1902850", # microtubule cytoskeleton organization involved in mitosis
                     "GO0007051", # spindle organization
                     "GO0051383", # kinetochore organization
                     "GO0007059" # chromosome segregation
) 
LateM_terms <- c("GO0000910", # cytokinesis
                 "GO0051302" # regulation of cell division
)
PhaseTransition_terms <- "GO0044770" #cell cycle phase transition



heatmap_annotation_row$S_phase <- apply(heatmap_annotation_row[,S_terms], 1, function(x){
  ifelse(sum(x) > 0, 1, 0)
})

heatmap_annotation_row$EarlyMidM_phase <- apply(heatmap_annotation_row[,EarlyMidM_terms], 1, function(x){
  ifelse(sum(x) > 0, 1, 0)
})

heatmap_annotation_row$LateM_phase <- apply(heatmap_annotation_row[,LateM_terms], 1, function(x){
  ifelse(sum(x) > 0, 1, 0)
})


heatmap_annotation_row <- as.data.frame(apply(heatmap_annotation_row, 2, as.factor))
rownames(heatmap_annotation_row) <- lookupGeneName(dat.filtered, genotypeCellType_diff_test_sig_genes)

heatmap_annotation_row <- heatmap_annotation_row[, c("S_phase", "EarlyMidM_phase", "LateM_phase", "GO0042063", "GO0007219", "GO0042136", "GO0010718")]

heatmap_data <- log10(exprs(dat.filtered)[genotypeCellType_diff_test_sig_genes,] + 1) # log transform
heatmap_data <- t(apply(heatmap_data, 1, function(x){ # Normalize gene expr by dividing by SD
  x/max(x)
}))

rownames(heatmap_data) <- lookupGeneName(dat.filtered,rownames(heatmap_data))
heatmap_data_neurons <- heatmap_data[,rownames(pData(dat.filtered)[pData(dat.filtered)$celltype == "neuron",])]
heatmap_data_neurons <- heatmap_data_neurons[,order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data_neurons)))))]

heatmap_data_glia <- heatmap_data[,rownames(pData(dat.filtered)[!pData(dat.filtered)$celltype == "neuron",])]
heatmap_data_glia <- heatmap_data_glia[,order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data_glia)))))]

heatmap_data <- cbind(heatmap_data_neurons, heatmap_data_glia)

heatmap_data <- heatmap_data[order.dendrogram(as.dendrogram(hclust(dist(heatmap_data)))),]

#heatmap_data <- heatmap_data[order.dendrogram(as.dendrogram(hclust(dist(heatmap_data)))),order.dendrogram(as.dendrogram(hclust(dist(t(heatmap_data)))))]

heatmap_annotation_col <- pData(dat.filtered)[,c("celltype","genotype")]

writeLines(paste0("\"", c("GO0007059", "GO0006260", "GO0045787", "GO0044770", "GO0006281", 
                          "GO0045786", "GO1902850", "GO0007051", "GO0051383", "GO0071103", 
                          "GO0051321", "GO0033044", "GO0051656", "GO0006310", "GO0034502", 
                          "GO0000910", "GO0051302", "GO0007219", "GO0042136", "GO0010718",
                          "GO0042063"), "\" = c(\"1\" = \"black\", \"0\" = \"grey90\")"), sep = ",\n")

heatmap_colors <- list(
  "Cell Type" =  c("neuron" = "#E69F00", "progenitor/glia" = "#56B4E9"),
  "Genotype" = c("het" = "#E41A1C", "hom" = "#377EB8"),
  
  "S phase" = c("1" = brewer.pal(7, "Dark2")[1], "0" = "grey90"),
  "early-mid M phase" = c("1" = brewer.pal(7, "Dark2")[2], "0" = "grey90"),
  "late M phase" = c("1" = brewer.pal(7, "Dark2")[3], "0" = "grey90"),
  "gliogenesis" = c("1" = brewer.pal(7, "Dark2")[4], "0" = "grey90"),
  "Notch signaling" = c("1" = brewer.pal(7, "Dark2")[5], "0" = "grey90"),
  "neurotransmitter biosynthetic process" = c("1" = brewer.pal(7, "Dark2")[6], "0" = "grey90"),
  "epithelial to mesenchymal transition" = c("1" = brewer.pal(7, "Dark2")[7], "0" = "grey90")
)

colnames(heatmap_annotation_col) <- c("Cell Type","Genotype")
colnames(heatmap_annotation_row) <- c("S phase", "early-mid M phase", "late M phase", "gliogenesis", "Notch signaling", "neurotransmitter biosynthetic process", "epithelial to mesenchymal transition")

pdf("SF6A-DEG_all_cells_genotypeXcelltype_heatmap.pdf", height = 20, width = 20)
pheatmap(mat = heatmap_data,
         scale="none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         drop_levels = FALSE,
         cluster_rows = F,
         cluster_cols = F,
         color = viridis(100),
         annotation_col = heatmap_annotation_col,
         annotation_row = heatmap_annotation_row,
         annotation_colors = heatmap_colors,
         annotation_names_col = TRUE,
         annotation_legend = TRUE
)
dev.off()


normalized_projections <- apply(normalized_projections, 1,
                                function(x){
                                  x=x-min(x)
                                  x=x/max(x)
                                })

pdf("SF6C-ecdf_all_patterns_facet_cluster.pdf", useDingbats = F)
merge(pData(dat.filtered)[,c("age", "genotype", "cluster")], normalized_projections, by = 0) %>%
  melt(., id.vars = c("Row.names", "age", "genotype", "cluster")) %>%
  mutate(condition = paste(age, genotype)) %>%
  ggplot(., aes(value, color = condition, lty = condition)) + 
  stat_ecdf(geom = "step") +
  facet_wrap(~cluster) +
  ggtitle("Cumulative density of normalized projection weight by cluster") +
  labs(x = "Pattern Weight", y = "Cumulative Density") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = c(0, 1.001)) +
  scale_color_manual("Age, Genotype", values = c("E12.5 het" = "#E41A1C", "E12.5 hom" = "#377EB8", "E14.5 het" = "#E41A1C", "E14.5 hom" = "#377EB8")) +
  scale_linetype_manual("Age, Genotype", values = c("E12.5 het" = "solid", "E12.5 hom" = "solid", "E14.5 het" = "dotted", "E14.5 hom" = "dotted")) +
  theme(aspect.ratio = 1, 
        strip.background = element_blank(),
        legend.position = "none", 
        plot.margin = margin(0, 2, 0, 0),
        axis.title = element_text(size = 8),
        axis.title.y = element_text(margin = margin(r = 0)),
        axis.title.x = element_text(margin = margin(t = -2)),
        axis.text = element_text(size = 6))

dev.off()

pdf("SF6D-ecdf_all_patterns.pdf", useDingbats = F, height = 3, width = 3)
merge(pData(dat.filtered)[,c("age", "genotype")], normalized_projections, by = 0) %>%
  melt(., id.vars = c("Row.names", "age", "genotype")) %>%
  mutate(condition = paste(age, genotype)) %>%
  ggplot(., aes(value, color = condition, lty = condition)) + 
  stat_ecdf(geom = "step") +
  ggtitle("All normalized projection weights") +
  labs(x = "Pattern Weight", y = "Cumulative Density") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = c(0, 1.001)) +
  scale_color_manual("Age, Genotype", values = c("E12.5 het" = "#E41A1C", "E12.5 hom" = "#377EB8", "E14.5 het" = "#E41A1C", "E14.5 hom" = "#377EB8")) +
  scale_linetype_manual("Age, Genotype", values = c("E12.5 het" = "solid", "E12.5 hom" = "solid", "E14.5 het" = "dotted", "E14.5 hom" = "dotted")) +
  theme(aspect.ratio = 1, 
        legend.position = c(0.5, 0.2),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 10),
        plot.margin = margin(0, 2, 0, 0),
        axis.title = element_text(size = 8),
        axis.title.y = element_text(margin = margin(r = 0)),
        axis.title.x = element_text(margin = margin(t = -2)),
        axis.text = element_text(size = 6))

dev.off()




#subset projection weights by condition (age*genotype)

tmp <- merge(pData(dat.filtered)[,c("cell_id", "age", "genotype", "cluster")] %>% 
               mutate(condition = paste(age, genotype)),
             normalized_projections, by.x = "cell_id", by.y = 0)

melted_tmp <- melt(tmp, id.vars = c("cell_id", "age", "genotype", "cluster", "condition"))

  


a <- melted_tmp[melted_tmp$condition == "E12.5 het", ]
b <- melted_tmp[melted_tmp$condition == "E12.5 hom", ]
c <- melted_tmp[melted_tmp$condition == "E14.5 het", ]
d <- melted_tmp[melted_tmp$condition == "E14.5 hom", ]

ks.test(a$value, b$value) # p = 0.0751
ks.test(a$value, c$value) # p < 2.2e-16 **
ks.test(a$value, d$value) # p = 1.665e-15 **
ks.test(b$value, c$value) # p < 2.2e-16 **
ks.test(b$value, d$value) # p = 1.388e-13 **
ks.test(c$value, d$value) # p = 1.315e-4 **

# 6 pairwise tests for each pattern
# * p < 0.05/6
# ** p < 0.01/6

calcKS <- function(cluster){
  if(length(a[a$cluster %in% cluster, "value"]) > 0 & length(b[b$cluster %in% cluster, "value"])){
    e <- ks.test(a[a$cluster %in% cluster, "value"], b[b$cluster %in% cluster, "value"]) # E12.5 het vs. E12.5 hom
    print(paste0("Cluster ", cluster, " E12.5 het vs. E12.5 hom: ", e$p.value))
  }
  if(length(a[a$cluster %in% cluster, "value"]) > 0 & length(c[c$cluster %in% cluster, "value"])){
    f <- ks.test(a[a$cluster %in% cluster, "value"], c[c$cluster %in% cluster, "value"]) # E12.5 het vs. E14.5 het
    print(paste0("Cluster ", cluster, " E12.5 het vs. E14.5 het: ", f$p.value))
  }
  if(length(a[a$cluster %in% cluster, "value"]) > 0 & length(d[d$cluster %in% cluster, "value"])){
    g <- ks.test(a[a$cluster %in% cluster, "value"], d[d$cluster %in% cluster, "value"]) # E12.5 het vs. E14.5 hom
    print(paste0("Cluster ", cluster, " E12.5 het vs. E14.5 hom: ", g$p.value))
  }
  if(length(b[b$cluster %in% cluster, "value"]) > 0 & length(c[c$cluster %in% cluster, "value"])){
    h <- ks.test(b[b$cluster %in% cluster, "value"], c[c$cluster %in% cluster, "value"]) # E12.5 hom vs. E14.5 het
    print(paste0("Cluster ", cluster, " E12.5 hom vs. E14.5 het: ", h$p.value))
  }
  if(length(b[b$cluster %in% cluster, "value"]) > 0 & length(d[d$cluster %in% cluster, "value"])){
    i <- ks.test(b[b$cluster %in% cluster, "value"], d[d$cluster %in% cluster, "value"]) # E12.5 hom vs. E14.5 hom
    print(paste0("Cluster ", cluster, " E12.5 hom vs. E14.5 hom: ", i$p.value))
  }
  if(length(c[c$cluster %in% cluster, "value"]) > 0 & length(d[d$cluster %in% cluster, "value"])){
    j <- ks.test(c[c$cluster %in% cluster, "value"], d[d$cluster %in% cluster, "value"]) # E14.5 het vs. E14.5 hom
    print(paste0("Cluster ", cluster, " E14.5 het vs. E14.5 hom: ", j$p.value))
  }
}

ks.test(a$Pattern_40, c$Pattern_40) # p < 2.2e-16 **, E12.5 het vs. E14.5 het (age)
ks.test(a$Pattern_40, d$Pattern_40) # p < 2.2e-16 **, E12.5 het vs E14.5 hom (age*genotype)
ks.test(b$Pattern_40, c$Pattern_40) # p < 2.2e-16 **, E12.5 hom vs. E14.5 het (age*genotype)
ks.test(b$Pattern_40, d$Pattern_40) # p = 3.331e-16 **, E12.5 hom vs. E14.5 hom (age)
ks.test(c$Pattern_40, d$Pattern_40) # p = 1.609e-5 **, E14.5 het vs. E14.5 hom (genotype)
