CC_patterns <- pData(dat.filtered)[,paste0("Linnarsson_Pattern_", c(47, 18, 27))]

CC_patterns <- apply(CC_patterns, 2, function(x){ #Normalize to [0,1]
  x=x-min(x)
  x=x/max(x)
})

CC_patterns <- melt(CC_patterns)
CC_patterns <- merge(pData(dat.filtered)[,c("cell_id", "UMAP1", "UMAP2")], CC_patterns, by.x = "cell_id", by.y = "Var1")
colnames(CC_patterns) <- c("cell_id", "x", "y", "pattern", "value")
CC_patterns$pattern <- str_remove(CC_patterns$pattern, "Linnarsson_")
CC_patterns$pattern <- factor(CC_patterns$pattern, levels = c("Pattern_47", "Pattern_18", "Pattern_27"))

phase <- c("Pattern_47" = "S", 
       "Pattern_18" = "Early M", 
       "Pattern_27" = "Late M")
#pdf("/Users/liz/Documents/Hopkins/Thesis Committee Meetings/AY 20-21/CC_pattern_projections.pdf",
#    height = 2.25,
#    width = 6,
#    useDingbats = F)
ggplot(CC_patterns) +
  geom_point(aes(x=x, y=y, 
                 color = value), size = 0.1) +
  facet_wrap(~pattern, scales = "free", nrow = 1, labeller = labeller(pattern = phase)) +
  scale_x_continuous("UMAP 1") +
  scale_y_continuous("UMAP2") +
  scale_color_viridis("Normalized\nProjection\nWeight", breaks = c(0,1)) +
  theme(aspect.ratio = 1,
        plot.margin = margin(0, 70, 0, 0),
        strip.text = element_text(margin = margin(b = 5)),
        strip.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(1, 0.4),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.height = unit(3,"mm"),
        legend.key.width = unit(3, "mm")
  ) 
#dev.off()


CC_phase <- pData(dat.filtered)[,paste0("Linnarsson_Pattern_", c(47, 18, 27))]

CC_phase <- as.data.frame(apply(CC_phase, 2, function(x){
  x=x-min(x)
  x=x/max(x)
}))

#m <- c("1" = "Early M",
#      "2" = "Late M")
  
#CC_phase$phase <- apply(CC_phase, 1, function(x){
#  ifelse(max(x) > 0.5,
#    ifelse(x[2] > 0.5 | x[3] > 0.5,
#      m[which.max(x[2:3])],
#      "S"),
#    "G1|G0")
#  })

m <- c("1" = "S",
  "2" = "Early M",
  "3" = "Late M")

CC_phase$phase <- apply(CC_phase, 1, function(x){
  ifelse(max(x) > 0.5,
         m[which.max(x)],
         "G1|G0")
  })

CC_phase <- merge(pData(dat.filtered)[,c("UMAP1", "UMAP2", "age", "sex", "genotype")], CC_phase, by = 0)

CC_phase$phase <- factor(CC_phase$phase, levels = c("G1|G0", "S", "Early M", "Late M"))

pdf("/Users/liz/Documents/Hopkins/Thesis Committee Meetings/AY 20-21/CC_phase.pdf", 
    height = 3.25, 
    width = 3.25,
    useDingbats = F)
ggplot(CC_phase) +
  geom_point(aes(x=UMAP1, y=UMAP2, 
                 color = phase), size = 0.5) +
  scale_color_manual(values = c("#969696", brewer.pal(3, "Set1"))) +
  guides(color = guide_legend(title = "CC Phase", override.aes = list(size = 3))) +
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.6, 0.22),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
  ) 
dev.off()

fd <- CC_phase
fd$phase <- factor(fd$phase, levels = c(levels(fd$phase), "S/G2/M"))
fd[fd$phase %in% c("S", "Early M", "Late M"), "phase"] <- "S/G2/M"

fd <- fd %>%
  group_by(age, genotype, phase) %>% 
  count() %>%
  group_by(age, genotype) %>%
  mutate(freq = n/sum(n)) %>% 
  mutate(condition = paste(genotype, phase))

#fd <- as.data.frame(fd)
#fd <- fd %>% 
#  add_row(age = "E14.5", genotype = "hom", phase = "S", n = 0, freq = 0, condition = "hom S") %>%
#  add_row(age = "E14.5", genotype = "hom", phase = "Early M", n = 0, freq = 0, condition = "hom Early M")

pdf("/Users/liz/Documents/Hopkins/Thesis Committee Meetings/AY 20-21/CC_proportion.pdf",
    height = 3.25,
    width = 4,
    useDingbats = F)
ggplot(fd, aes(x = age, y = freq, group = condition, color = phase, lty = genotype)) +
  geom_line() +
  geom_point() +
  scale_color_manual("CC Phase", values = c("black", "red")) +
  scale_y_continuous("Proportion", limits = c(0, 1), expand = c(0,0)) +
  scale_x_discrete("Age", expand = c(0.05, 0.05)) +
  scale_linetype("Genotype")
dev.off()

my_color_scale <- scale_color_gradientn("", colors = c("#F1ECB9", "#BEDEE5", "#25C6FE", "#9232B4", "#2C0404"), breaks = c(0, 1), labels = c("min", "max"))
my_theme <-   theme(strip.background = element_blank(),
                    strip.text = element_text(size = 20, margin = margin(b = 5)),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    #axis.line = element_blank(),
                    axis.title = element_blank())

pdf("/Users/liz/Documents/Hopkins/Thesis Committee Meetings/AY 20-21/RetCFP.pdf",
    height = 3.25,
    width = 5.75,
    useDingbats = F)
plotUMAP(dat.filtered, markers = c("Ret", "tCFP"), scale = T, size = 0.1) +
  my_color_scale +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 20, margin = margin(b = 5)),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.85, 0.3))
dev.off()

pdf("/Users/liz/Documents/Hopkins/Thesis Committee Meetings/AY 20-21/CFPCartptGal.pdf",
    height = 3.25,
    width = 8.25,
    useDingbats = F)
plotUMAP(dat.filtered, markers = c("tCFP", "Cartpt", "Gal"), scale = T, size = 0.1) +
  my_color_scale +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 20, margin = margin(b = 5)),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.9, 0.3))
dev.off()


pdf("/Users/liz/Documents/Hopkins/Thesis Committee Meetings/AY 20-21/CFPRetGfra1.pdf",
    height = 6,
    width = 6,
    useDingbats = F)
plotUMAP(dat.filtered, markers = c("tCFP", "Ret", "Gfra1", "Gfra2"), scale = T, size = 0.1) +
  my_color_scale +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 20, margin = margin(b = 5)),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.85, 0.15))
dev.off()

pdf("/Users/liz/Documents/Hopkins/Thesis Committee Meetings/AY 20-21/Plp1.pdf",
    height = 3.25,
    width = 3.25,
    useDingbats = F)
plotUMAP(dat.filtered, markers = c("Plp1"), scale = T, size = 0.1, nrow = 1) +
  my_color_scale +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 20, margin = margin(b = 5)),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.75, 0.3))
dev.off()
