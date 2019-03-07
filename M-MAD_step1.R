#!/usr/bin/env Rscript
# script to perform M-MAD analysis 

##' @param data data.frame of G-MAD results containing -log10(p values) with pathways as rows and genes as columns, obtained from "peer.R". 1st row contains pathway ids. 
##' @param pathways list of pathways containing the gene entrez id for each pathway
##' @return results.all data.table of camera p values between pathways and pathways. 


mmad_step1 <- function(data, pathways){
  require(limma)
  require(data.table)
  options(stringsAsFactors = F)
  set.seed(666)
  
  # check for missing values
  if(sum(is.na(data)) != 0){
    message('there are NA in the data')
    # remove rows containing NAs
    data <- data[which(rowSums(is.na(data)) == 0), ]
  }
  
  # create common index
  index <- ids2indices(pathways, colnames(data)[-1], remove.empty=FALSE)
  
  for(i in 1:nrow(data)){
    pathway <- data[i, 'pathway']
    print(c(i, pathway))
    data.i <- as.numeric(data[i, -1])
    
    stats <- setNames(data.i, colnames(data)[-1])
    result.i <- tryCatch(cameraPR(statistic = stats, index = index, inter.gene.cor=0.01, use.ranks = F, sort = F), error=function(e) FALSE)
    if(is.logical(result.i)) next
    
    result.i$path_id <- rownames(result.i)
    result.i$logp <- round(-log10(result.i$PValue), digits = 3)
    result.i[which(result.i$Direction == "Down"),"logp"] <- -result.i[which(result.i$Direction == "Down"),"logp"]
    result.i <- result.i[,which(colnames(result.i) %in% c("path_id", "logp"))]
    colnames(result.i)[which(colnames(result.i) %in% c("logp"))] <- pathway
    
    result.i <- data.table(result.i, key = "path_id")
    if(i == 1){
      result.all <- result.i
    }else{
      result.all <- merge(result.all, result.i, all = TRUE) 
    }
  }
  
  return(result.all)
}

