#!/usr/bin/env Rscript
# script to run pathway analysis for all genes using camera 

##' @param arrayData data.frame of expression residual with genes as rows and individuals as columns, obtained from "peer.R"
##' @param pathways list of pathways containing the gene entrez id for each pathway
##' @return results.all data.table of camera p values between genes and pathways with pathways as rows and genes as columns

gmad_step1 <- function(arrayData, pathways){
  require(limma)
  require(data.table)
  options(stringsAsFactors = F)
  set.seed(666)
  
  arrayData <- as.matrix(arrayData)
  
  # generate the index vector, made using ids2indices
  index <- ids2indices(pathways, rownames(arrayData), remove.empty=F)
  
  all.gene.num <- nrow(arrayData)
  for(target.gene.i in 1:all.gene.num){
    target.gene <- rownames(arrayData)[target.gene.i]
    message(paste0(target.gene.i, '/',all.gene.num))
    gene.i <- as.numeric(arrayData[target.gene.i, ])
    
    if(sd(gene.i) == 0) { # if gene.i has unique values, proceed with next gene
      message(paste(target.gene, "contains unique values, proceed with next gene")) 
      next
    }
    
    design <- model.matrix( ~ gene.i)
    results <- camera(arrayData, index, design, inter.gene.cor=0.01, allow.neg.cor = T, use.ranks = F, sort = F)
    
    # !!! here two-sided p-values are saved, should be converted to one-sided p-value (all positive side) for calculation
    results$PValue <- results$PValue / 2        # convert to one-sided pvalues (divide 2-side p by 2)
    
    results$pathway <- rownames(results)
    results$logp <- round(-log10(results$PValue), digits = 3)
    results[which(results$Direction == "Down"),"logp"] <- -results[which(results$Direction == "Down"),"logp"]
    results.x <- results[,which(colnames(results) %in% c("pathway", "logp"))]
    colnames(results.x)[which(colnames(results.x) %in% c("logp"))] <- target.gene
    results.x <- data.table(results.x, key = "pathway")
    
    if(target.gene.i == 1){
      results.all <- results.x
    }else{
      results.all <- results.all[results.x]
    }
  }
  return(results.all)
}


