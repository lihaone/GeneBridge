#!/usr/bin/env Rscript
# compute the correlation coefficient between genes in respective pathways

##' @param arrayData data.frame of expression residual with genes as rows and individuals as columns, obtained from "peer.R"
##' @param pathways list of pathways containing the gene entrez id for each pathway, obtained from "load.pathways" function in "utils.R"
##' @return r.mean data.frame of the average correlation coefficients of all pathways

average_correlation_coefficients <- function(arrayData, pathways){
  rmatrix <- cor(t(arrayData), use='pairwise.complete.obs', method='pearson')
  rmatrix <- round(rmatrix, digits=3)

  r.mean <- data.frame(pathway=names(pathways), gene.num=NA, mean=NA)
  for(i in 1:length(pathways)){
    pathway <- names(pathways)[i]
    print(c(i, pathway))
    match.genes <- which(rownames(rmatrix) %in% pathway)
    r.mean[i, 'gene.num'] <- length(match.genes)
    r.mean[i, 'mean'] <- mean(rmatrix[match.genes, match.genes], na.rm=T)
  }
  return(r.mean)
}
