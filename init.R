library(reticulate)
use_python(python="/Library/Frameworks/Python.framework/Versions/2.7/bin/python",required=TRUE)
library(umapr)
library(ggplot2)
library(RColorBrewer)
library(ggbiplot)
library(monocle)
library(gridExtra)
library(reshape2)
library(tidyverse)
library(devtools)
library(Matrix)
library(CoGAPS)
library(stringr)
library(cowplot)
library(viridis)
library(corrplot)
library(projectR)

#add shape arg to UMAP
plotUMAP <- function(cds, markers = NULL, color_by = NULL, shape = NULL, logMode = T, scaled = F, size = 1.5, ...){
  pheno <- pData(cds)
  expr <- exprs(cds)
  features <- fData(cds)
  if(!is.null(markers)){
    if(sum(!markers %in% features$gene_short_name) > 0){
      warning(paste(paste0(markers[!markers %in% features$gene_short_name], collapse = ", " ), "is/are not valid gene names"))
      markers <- markers[markers %in% features$gene_short_name]
    }
    genes <- expr[rownames(features) %in% lookupGeneId(cds, markers), , drop = F]
    if(logMode){
      genes <- log10(genes + 1)
    }
    if(scaled){
        geneMax <- rowMax(genes)
        genes <- genes/geneMax
    }
    genes <- t(genes)
    genes <- melt(genes)
    colnames(genes) <- c("cell_id", "gene_id", "value")
    genes <- merge(genes, features[, c("gene_id", "gene_short_name")], by.x = "gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
    tmp <- merge(pheno, genes, by.x = 0, by.y = "cell_id", sort = FALSE)
    if(!is.null(shape)){
      shape <- rep(pheno[, shape], length(markers))
    } 
    p <- ggplot(tmp, aes(x = UMAP1, y = UMAP2)) + 
      geom_point(aes_string(shape = shape, color = "value"), size = size, ...)
    if(logMode){
      if(scaled){
        p <- p + scale_color_viridis("log10(CPC + 1)\nscaled to max", option = "C")   
      }else{
        p <- p + scale_color_viridis("log10(CPC + 1)", option = "C")
      }
    }else{
      if(scaled){
        p <- p + scale_color_viridis("CPC\nscaled to\nmax", option = "C")
      }else{
      p <- p + scale_color_viridis("CPC", option = "C")
      }
    }
    if(length(markers) > 1){
      p <- p + facet_wrap('gene_short_name')
    }else{
      p <- p + ggtitle(markers)
    }
  }else{
    p <- ggplot(pheno, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes_string(color = color_by, shape = shape), size = size, ...)
    if(is.discrete(color_by)){
      n_colors <- length(unique(pheno[,color_by]))
      if(n_colors <= 9){
        p <- p + scale_color_brewer(palette = "Set1")
      }else{
        p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(n_colors))
      }
    }else{
      p <- p + scale_color_viridis(option = "C")
    }
  }
  return(p)
}

##################
#Helper Utils
####################
lookupGeneId<-function(eset,gene_names, unique = T){
  res <- rownames(fData(eset))[fData(eset)$gene_short_name %in% gene_names]
  if(unique){
  res <- unique(res)
  }
  res
}

lookupGeneName<-function(eset, gene_id, unique = T){
  gene_id <- gsub("\\..*", "", gene_id)
  rownames(fData(eset)) <- gsub("\\..*", "", rownames(fData(eset)))
  res <- fData(eset)[gene_id, "gene_short_name"]
  if(unique){
    res <- unique(res)
  }
  res
}

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
