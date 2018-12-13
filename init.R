library(reticulate)
use_python(python="/Library/Frameworks/Python.framework/Versions/2.7/bin/python",required=TRUE)
library(umapr)
library(ggplot2)
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

plotUMAP <- function(cds, markers = NULL, logMode = T, color_by = NULL, scaled = F, sort = F, size = 1.5){
  pheno <- pData(cds)
  expr <- exprs(cds)
  features <- fData(cds)
  if(!is.null(markers)){
    genes <- expr[rownames(features) %in% lookupGeneId(cds, markers), , drop = F]
    if(scaled){
        geneMeans <- rowMeans(genes)
        #genes <- apply(genes, 1, function(x){(x - mean(x))/sd(x)})
        genes <- genes/geneMeans
    }
    genes <- t(genes)
    if(logMode){
      genes <- log10(genes + 1)
    }
    genes <- melt(genes)
    colnames(genes) <- c("cell_id", "gene_id", "value")
    genes <- merge(genes, features, by.x = "gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
    tmp <- merge(pheno, genes, by.x = 0, by.y = "cell_id", sort = FALSE)
    p <- ggplot(tmp, aes(x = UMAP1, y = UMAP2)) + 
      geom_point(aes_string(color = "value"), size = size) + 
      theme_bw()
    if(logMode){
      if(scaled){
        p <- p + scale_color_viridis("log10(value + 1)\nscaled to mean", option = "C")   
      }else{
        p <- p + scale_color_viridis("log10(value + 1)", option = "C")
      }
    }else{
      if(scaled){
        p <- p + scale_color_viridis("scaled to mean", option = "C")
      }else{
      p <- p + scale_color_viridis(option = "C")
      }
    }
    if(length(markers) > 1){
      p <- p + facet_wrap('gene_short_name')
    }else{
      p <- p + ggtitle(markers)
    }
  }else{
    p <- ggplot(pheno, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes_string(color = color_by), size = size) + 
      theme_bw()
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


# fix patternMarker() bug
patternMarkers <- function (Amatrix = NA, scaledPmatrix = FALSE, Pmatrix = NA, 
          threshold = "all", lp = NA, full = FALSE) 
{
  if (scaledPmatrix == FALSE) {
    if (!is.na(Pmatrix)) {
      pscale <- apply(Pmatrix, 1, max)
      Amatrix <- sweep(Amatrix, 2, pscale, FUN = "*")
    }
    else {
      warning("P values must be provided if not already scaled")
    }
  }
  Arowmax <- t(apply(Amatrix, 1, function(x) x/max(x)))
  if (!is.na(lp)) {
    if (length(lp) != dim(Amatrix)[2]) {
      warning("lp length must equal the number of columns of the Amatrix")
    }
    sstat <- apply(Arowmax, 1, function(x) sqrt(t(x - lp) %*% 
                                                  (x - lp)))
    ssranks <- rank(sstat)
    ssgenes.th <- names(sort(sstat, decreasing = FALSE, na.last = TRUE))
  }
  else {
    sstat <- matrix(NA, nrow = nrow(Amatrix), ncol = ncol(Amatrix), 
                    dimnames = dimnames(Amatrix))
    ssranks <- matrix(NA, nrow = nrow(Amatrix), ncol = ncol(Amatrix), 
                      dimnames = dimnames(Amatrix))
    ssgenes <- matrix(NA, nrow = nrow(Amatrix), ncol = ncol(Amatrix), 
                      dimnames = NULL)
    nP<-dim(Amatrix)[2]
    for (i in 1:nP) {
      lp <- rep(0, dim(Amatrix)[2])
      lp[i] <- 1
      sstat[, i] <- unlist(apply(Arowmax, 1, function(x) sqrt(t(x - 
                                                                  lp) %*% (x - lp))))
      ssranks[, i] <- rank(sstat[, i])
      ssgenes[, i] <- names(sort(sstat[, i], decreasing = FALSE, 
                                 na.last = TRUE))
      ssgenes[, i] <- gsub("\\..*","", ssgenes[, i])
    }
    if (threshold == "cut") {
      geneThresh <- sapply(1:nP, function(x) min(which(ssranks[ssgenes[, 
                                                                       x], x] > apply(ssranks[ssgenes[, x], ], 1, min))))
      ssgenes.th <- sapply(1:nP, function(x) ssgenes[1:geneThresh[x], 
                                                     x])
    }
    else if (threshold == "all") {
      pIndx <- unlist(apply(sstat, 1, which.min))
      gBYp <- list()
      for (i in sort(unique(pIndx))) {
        gBYp[[i]] <- sapply(strsplit(names(pIndx[pIndx == 
                                                   i]), "[.]"), function(x) x[[1]][1])
      }
      ssgenes.th <- lapply(1:max(sort(unique(pIndx))), 
                           function(x) ssgenes[which(ssgenes[, x] %in% gBYp[[x]]), 
                                               x])
    }
    else {
      stop("Threshold arguement not viable option")
    }
  }
  if (full) {
    return(list(PatternMarkers = ssgenes.th, PatternRanks = ssranks, 
                PatternMarkerScores = sstat))
  }
  else {
    return(PatternMarkers = ssgenes.th)
  }
}

# fix calcCoGAPSStat() bug

calcCoGAPSStat <- function (Amean, Asd, GStoGenes, numPerm = 500) 
{
  if (sum(Asd == 0) > 0 | sum(is.na(Asd)) > 0) 
    Asd[Asd == 0 | is.na(Asd)] <- 1e-06
  zMatrix <- calcZ(Amean, Asd)
  if (!is(GStoGenes, "data.frame") && !is(GStoGenes, "list") && 
      !is(GStoGenes, "GSA.genesets")) {
    stop("GStoGenes must be a data.frame,GSA.genesets, or list with format specified in the users manual.")
  }
  if (is(GStoGenes, "GSA.genesets")) {
    names(GStoGenes$genesets) <- GStoGenes$geneset.names
    GStoGenes <- GStoGenes$genesets
  }
  if (is(GStoGenes, "list")) {
    GStoGenesList <- GStoGenes
  }
  else {
    GStoGenesList <- list()
    for (i in 1:dim(GStoGenes)[2]) {
      GStoGenesList[[as.character(colnames(GStoGenes)[i])]] <- as.character(unique(GStoGenes[, 
                                                                                             i]))
    }
  }
  numGS <- length(names(GStoGenesList))
  numPatt <- dim(zMatrix)[2]
  numG <- dim(zMatrix)[1] + 0.9999
  statsUp <- matrix(ncol = numGS, nrow = numPatt)
  statsDown <- matrix(ncol = numGS, nrow = numPatt)
  actEst <- matrix(ncol = numGS, nrow = numPatt)
  results <- list()
  zPerm <- matrix(ncol = numPerm, nrow = numPatt)
  for (gs in 1:numGS) {
    genes <- GStoGenesList[[names(GStoGenesList)[gs]]]
    index <- gsub("\\..*","", rownames(zMatrix)) %in% genes
    zValues <- zMatrix[index, 1]
    numGenes <- length(zValues)
    label <- as.character(numGenes)
    if (!any(names(results) == label)) {
      for (p in 1:numPatt) {
        for (j in 1:numPerm) {
          temp <- floor(runif(numGenes, 1, numG))
          temp2 <- zMatrix[temp, p]
          zPerm[p, j] <- mean(temp2)
        }
      }
      results[[label]] <- zPerm
    }
  }
  for (p in 1:numPatt) {
    for (gs in 1:numGS) {
      genes <- GStoGenesList[[names(GStoGenesList)[gs]]]
      index <- gsub("\\..*","", rownames(zMatrix)) %in% genes
      zValues <- zMatrix[index, p]
      zScore <- mean(zValues)
      numGenes <- length(zValues)
      label <- as.character(numGenes)
      permzValues <- results[[label]][p, ]
      ordering <- order(permzValues)
      permzValues <- permzValues[ordering]
      statistic <- sum(zScore > permzValues)
      statUpReg <- 1 - statistic/length(permzValues)
      statsUp[p, gs] <- max(statUpReg, 1/numPerm)
      statistic <- sum(zScore < permzValues)
      statDownReg <- 1 - statistic/length(permzValues)
      statsDown[p, gs] <- max(statDownReg, 1/numPerm)
      activity <- 1 - 2 * max(statUpReg, 1/numPerm)
      actEst[p, gs] <- activity
    }
  }
  colnames(statsUp) <- names(GStoGenesList)
  colnames(statsDown) <- names(GStoGenesList)
  colnames(actEst) <- names(GStoGenesList)
  rownames(statsUp) <- colnames(zMatrix)
  rownames(statsDown) <- colnames(zMatrix)
  rownames(actEst) <- colnames(zMatrix)
  results[["GSUpreg"]] <- statsUp
  results[["GSDownreg"]] <- statsDown
  results[["GSActEst"]] <- actEst
  return(results)
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

lookupGeneName<-function(eset,gene_id, unique = T){
  res <- fData(eset[gene_id,])$gene_short_name
  if(unique){
    res <- unique(res)
  }
  res
}