gene_ids <- str_split(ego_glia_genotype@result[glia_terms, "geneID"], "\\/")

names(gene_ids) <- ego_glia_genotype@result[glia_terms, "ID"]

redundant <- list()
for(i in seq_along(gene_ids)){
#for(i in 1:5){
  for(j in seq_along(gene_ids)){
#  for(j in 1:5){
    if(i != j & sum(gene_ids[[i]] %in% gene_ids[[j]])/length(gene_ids[[i]]) > 0.9){
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
ego_genotype_subset <- ego_genotype_subset[!ego_genotype_subset$ID %in% redundant,]

gene_ids <- str_split(ego_genotype_subset$geneID, "\\/")

names(gene_ids) <- ego_genotype_subset$ID

redundant <- list()
for(i in seq_along(gene_ids)){
  #for(i in 1:5){
  for(j in seq_along(gene_ids)){
    #  for(j in 1:5){
    if(i != j & sum(gene_ids[[i]] %in% gene_ids[[j]])/length(gene_ids[[i]]) > 0.80){
      redundant[length(redundant) + 1] <- names(gene_ids)[i]
      break
#            print(paste("i =", i, "j =", j))
#            print(gene_ids[[i]] %in% gene_ids[[j]])
#            print(paste("names(gene_ids)[[i]] =", names(gene_ids)[[i]]))
#            print(sum(gene_ids[[i]] %in% gene_ids[[j]])/length(gene_ids[[i]]))
#            break
    }
  }
}

redundant <- unlist(redundant)
ego_genotype_subset <- ego_genotype_subset[!ego_genotype_subset$ID %in% redundant,]
