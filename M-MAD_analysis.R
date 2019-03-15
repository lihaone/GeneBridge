#!/usr/bin/env Rscript
# script to perform M-MAD analysis 


# example:
species <- 'human'
data.files <- list.files('./data/output/GMAD_module_preBonf/', full.names=T)
pathways <- load.pathways(dir="./data/utils data/pathway/", species='human')
out.dir <- './data/output/MMAD_module/'

data.file <- data.files[1]
print(data.file)
pathway <- sapply(strsplit(data.file, split='//', fixed=TRUE), function(x) (x[2]))
pathway <- sapply(strsplit(pathway, split='--module_', fixed=TRUE), function(x) (x[2]))
pathway <- gsub(pattern=".gz", replacement="", x=pathway, fixed=T)

data <- read.table(data.file, header=T, sep='\t', quote="", row.names=1)

result.all <- mmad_step1(data, pathways)
result.all.sig <- mmad_step2(result.all, pathways, sample.size, r_mean)

write.table(result.all.sig, gzfile(paste0(out.dir,species,'--module_',pathway,'.gz')), quote=F, sep="\t", row.names=F, col.names=T)

##' @param data G-MAD results obtained from "gmad_step2_module" function in "G-MAD_step2.R". loaded from files in folder "./data/output/GMAD_module_preBonf/"
##'             data.frame containing -log10(p values) with pathways as rows and genes as columns, 1st row contains pathway ids. 
##' @param pathways list of pathways containing the gene entrez id for each pathway, obtained from "load.pathways" function in "utils.R"
##' @return result.all data.table of camera p values between pathways and pathways. 
mmad_step1 <- function(data, pathways){
  require(limma)
  require(data.table)
  set.seed(666)
  
  result.all <- NULL
  for(i in 1:ncol(data)){
    dataset.id <- colnames(data)[i]
    print(c(i, dataset.id))
    data.i <- data.frame(gene_id=rownames(data), value=as.numeric(data[, i]))
    data.i <- data.i[complete.cases(data.i), ]
    
    # create index
    index <- ids2indices(pathways, data.i$gene_id, remove.empty=FALSE)
    
    stats <- setNames(data.i$value, data.i$gene_id)
    result.i <- tryCatch(cameraPR(statistic=stats, index=index, inter.gene.cor=0.01, use.ranks=F, sort=F), error=function(e) FALSE)
    if(is.logical(result.i)) next
    
    result.i$path_id <- rownames(result.i)
    result.i$logp <- round(-log10(result.i$PValue), digits=3)
    result.i[which(result.i$Direction == "Down"),"logp"] <- -result.i[which(result.i$Direction == "Down"),"logp"]
    result.i <- result.i[,which(colnames(result.i) %in% c("path_id", "logp"))]
    result.i$logp[is.infinite(result.i$logp)] <- max(result.i$logp[is.finite(result.i$logp)], na.rm=T)  # replace Inf from -log10(p) by the maximum of the rest values
    colnames(result.i)[which(colnames(result.i) %in% c("logp"))] <- dataset.id
    
    result.i <- data.table(result.i, key="path_id")
    if(is.null(result.all)){
      result.all <- result.i
    }else{
      result.all <- merge(result.all, result.i, all=TRUE) 
    }
  }
  return(result.all)
}


##' @param result.all results obtained from "mmad_step1" function in "M-MAD_analysis.R". 
##' @param pathways list of pathways/modules containing the gene entrez id for each pathway/module, obtained from "load.pathways" function in "utils.R"
##' @param sample.size data.frame containing the sample size information for all datasets, obtained from "load.sample.size" function in "utils.R"
##' @param r_mean data.frame containing the average correlation coefficients for pathways/modules in all datasets, obtained from "load.r_mean" function in "utils.R"
##' @return result.all.sig data.frame of final association between pathways/modules and pathways/modules. 
mmad_step2 <- function(result.all, pathways, sample.size, r_mean){
  result.all.sig <- result.all
  class(result.all.sig) <- "data.frame"
  rownames(result.all.sig) <- result.all.sig[,1]
  result.all.sig <- result.all.sig[, -1]
  
  logp.thres <- -log10(0.05/nrow(result.all.sig))
  result.all.sig <- data.matrix(result.all.sig)
  result.all.sig[which((result.all.sig > -logp.thres) & (result.all.sig < logp.thres))] <- 0
  result.all.sig[which(result.all.sig >= logp.thres)] <- 1
  result.all.sig[which(result.all.sig <= -logp.thres)] <- -1
  result.all.sig <- data.frame(result.all.sig)
  
  ### remove datasets for pathways based on enrichment results against themselves
  data.row <- which(rownames(result.all.sig) == pathway)
  if(length(data.row) != 0){
    cols.rm <- which(result.all.sig[data.row, ] == 0)
    if(length(cols.rm) != 0){
      print(paste0('removing ', length(cols.rm), ' columns for ', pathway))
      result.all.sig <- result.all.sig[, -cols.rm]
    }
  }
  
  colnames(result.all.sig) <- make.names(colnames(result.all.sig))
  sample.size$tissue <- make.names(sample.size$tissue)
  colnames(r_mean) <- make.names(colnames(r_mean))
  datasets <- intersect(intersect(sample.size$tissue, colnames(result.all.sig)), colnames(r_mean)) # get the datasets with all data
  
  # compute the weighted average
  sample.size.use <- sample.size[which(sample.size$tissue %in% datasets),] # should only use the common ones in match
  
  # get the gene correlations for the pathway and matched data.files
  r_mean.select <- r_mean[which(rownames(r_mean) == pathway), which(colnames(r_mean) %in% datasets)]
  r_mean.use <- data.frame(tissue=colnames(r_mean.select), r.mean=as.numeric(r_mean.select[1,]))
  sample.size.use <- plyr::join(sample.size.use, r_mean.use, by='tissue', type='inner')
  sample.size.use <- sample.size.use[complete.cases(sample.size.use), ]
  result.all.sig <- result.all.sig[, match(sample.size.use$tissue, colnames(result.all.sig))]
  
  if(is.null(ncol(result.all.sig))) next
  if(ncol(result.all.sig) == 0) next
  
  result.all.sig$mean.w <- apply(result.all.sig, 1, FUN=function(x) weighted.mean(x, w=sqrt(sample.size.use$size)*sample.size.use$r.mean, na.rm=T))
  result.all.sig$mean.w <- round(result.all.sig$mean.w, digits=3)
  result.all.sig$path_id <- rownames(result.all.sig)
  result.all.sig <- result.all.sig[,c(ncol(result.all.sig), 1:(ncol(result.all.sig)-1))]
  
  return(result.all.sig)
}
