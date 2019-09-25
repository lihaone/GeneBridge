#!/usr/bin/env Rscript
# script to remove both known and hidden covariates using PEER
# partially adapted from: https://github.com/hall-lab/gtex

##' @param mat numeric data.matrix of gene expressions with genes as rows and individuals as columns, with gene entrez ids as rownames
##' @param mat.aligner data.frame of 2 columns contaning the gene entrez id and gene symbols, with "id" and "gene" as colnames
##' @param covs data.frame contaning the meta data / covariate data,
##' @return residuals data.frame of expression residual with genes as rows and individuals as columns


peer.residual <- function(mat, mat.aligner, covs=NULL, n.iteration=1000){
  library(peer)
  set.seed(666)

  model <- PEER()

  # normalize the values of one gene to the range of 0~1, so that no variable's range has an unduly large impact on the distance measurement.
  # ref: https://stackoverflow.com/questions/16874038/error-with-knn-function
  range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  mat <- t(apply(mat, 1, FUN=function(x) range01(x, na.rm=T)))

  # optional: Impute missing values in gene expression
  if(sum(is.na(mat)) != 0){
    # remove genes with more than 30% missing values
    row.NAs <- which(rowSums(is.na(mat)) >= 0.3 * ncol(mat))
    if(length(row.NAs) != 0){
      message("removing genes with more than 30% missing values")
      mat <- mat[-row.NAs, ]
    }

    # remove samples with more than 30% missing values
    col.NAs <- which(colSums(is.na(mat)) >= 0.3 * nrow(mat))
    if(length(col.NAs) != 0){
      message("removing samples with more than 30% missing values")
      mat <- mat[, -col.NAs]
    }

    if(sum(is.na(mat)) != 0){
      message("imputing missing values in gene expression data")
      library(impute)
      mat.impute <- impute.knn(mat)
      mat <- mat.impute$data
    }
  }

  PEER_setPhenoMean(model, t(as.matrix(mat)))

  # set covariates
  if(!is.null(covs)){
    PEER_setCovariates(model, as.matrix(covs))
  }

  # set the number of covariates 25% of the sample number, but less than 100. (just a sufficiently large value)
  # ref: (PEER paper) https://www.nature.com/articles/nprot.2011.457
  numcov <- 0.25*ncol(mat)
  numcov <- min(numcov, 100)

  PEER_setNk(model, numcov)
  PEER_getNk(model)

  # use the default tolerances
  PEER_setTolerance(model, 0.001)
  PEER_setVarTolerance(model, 0.00001)

  # set the number of iterations
  PEER_setNmax_iterations(model, n.iteration)

  PEER_update(model)

  factors=t(PEER_getX(model))
  weights=PEER_getW(model)
  precision=PEER_getAlpha(model)

  residuals=t(PEER_getResiduals(model))
  rownames(residuals) <- rownames(mat)
  colnames(residuals) <- colnames(mat)
  residuals <- as.data.frame(round(residuals, digits=4))

  residuals <- cbind(rownames(residuals), residuals)
  colnames(residuals)[1] <- "id"
  residuals <- plyr::join(residuals, mat.aligner, by="id")
  residuals <- residuals[,c(1,ncol(residuals),2:(ncol(residuals)-1))]

  return(residuals)
}
