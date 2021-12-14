gene_ids <- str_split(ego_glia_genotype@result[ego_glia_genotype@result$qvalue < 0.05, "geneID"], "\\/")

names(gene_ids) <- ego_glia_genotype@result[ego_glia_genotype@result$qvalue < 0.05, "ID"]

redundant <- list()
#for(i in seq_along(gene_ids)){
for(i in 1:10){
  for(j in seq_along(gene_ids)){
#  for(j in 1:5){
    if(i != j & sum(gene_ids[[i]] %in% gene_ids[[j]])/length(gene_ids[[i]]) == 1){
      redundant[length(redundant) + 1] <- names(gene_ids)[i]
            print(paste("i =", i, "j =", j))
            print(gene_ids[[i]] %in% gene_ids[[j]])
            print(paste("names(gene_ids)[[i]] =", names(gene_ids)[[i]]))
            print(sum(gene_ids[[i]] %in% gene_ids[[j]])/length(gene_ids[[i]]))
      break

    }
  }
}

redundant <- unlist(redundant)
ego_glia_genotype_subset1 <- ego_glia_genotype@result[ego_glia_genotype@result$qvalue < 0.05 & !ego_glia_genotype@result$ID %in% redundant,]
