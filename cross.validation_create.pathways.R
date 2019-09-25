# create query pathways for cross validation
### load pathways
set.seed(666)
pathway.dir <- './data/utils data/pathway/'
all.species <- c('human', 'mouse', 'rat', 'fly', 'worm', 'yeast')

for(species.i in all.species){
  # species.i <- 'human'
  pathways <- load.pathways(dir=pathway.dir, species=species.i)

  # output variables
  pathways.new <- list()
  path.genes.left <- list()

  for(path.i in 1:length(pathways)){
    pathway <- names(pathways)[path.i]
    print(c(path.i, pathway))
    pathway.genes <- pathways[[path.i]]

    if(length(pathway.genes) <= 50){           # use leave-one-out cross-validation
      for(i in 1:length(pathway.genes)){
        gene.i <- pathway.genes[i]
        # print(c(i, gene.i))

        pathway.new.i <- pathway.genes[-i]
        pathway.new.i <- list(pathway.new.i)
        path.genes.left.i <- list(pathway.genes[i])

        names(pathway.new.i) <- names(path.genes.left.i) <- paste0(pathway,'__LOOCV_',i)

        pathways.new <- c(pathways.new, pathway.new.i)
        path.genes.left <- c(path.genes.left, path.genes.left.i)
      }
    }else{                                     # use 10-fold cross-validation
      library(caret)
      pathway.new.i <- pathway.index.i <- createFolds(pathway.genes, k=10, list=T, returnTrain=T)
      if(length(pathway.new.i) != 10){
        pathway.new.i <- pathway.index.i <- createFolds(pathway.genes, k=10, list=T, returnTrain=T)
      }
      path.genes.left.i <- pathway.new.i

      for(i in 1:length(pathway.index.i)){
        # print(i)
        pathway.new.i[[i]] <- pathway.genes[pathway.index.i[[i]]]
        path.genes.left.i[[i]] <- pathway.genes[-pathway.index.i[[i]]]
      }
      names(path.genes.left.i) <- names(pathway.new.i) <- paste0(pathway,'__10fold_CV_',1:10)
      pathways.new <- c(pathways.new, pathway.new.i)
      path.genes.left <- c(path.genes.left, path.genes.left.i)
    }
    print(paste0('pathway numbers: ', length(pathways.new)))
    print(paste0('left gene numbers: ', length(path.genes.left)))
  }

  path.genes.left.df <- melt(path.genes.left)
  colnames(path.genes.left.df) <- c('gene_id', 'pathway')

  saveRDS(pathways.new, paste0(pathway.dir, species.i, '_cross_validation_all_pathway.RDS'))
  saveRDS(path.genes.left.df, paste0(pathway.dir, species.i, '_cross_validation_left_pathway_genes.RDS'))
}
