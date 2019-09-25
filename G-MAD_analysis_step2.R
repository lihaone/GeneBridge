#!/usr/bin/env Rscript
# second step of G-MAD
# meta-analysis of results obtained from 'G-MAD_analysis_step1.R' for each gene across all analyzed datasets

################################################################################################################################################
######  analysis for one gene vs all modules/pathways
##' @param data.files data files in RDS format for all datasets resulted from 'G-MAD_analysis_step1.R'
##' @param sample.size data.frame containing the sample size information for all datasets, obtained from 'load.sample.size' function in 'utils.R'
##' @param r_mean data.frame containing the average correlation coefficients for pathways in all datasets, obtained from 'load.r_mean' function in 'utils.R'
##' @param gene2pathway list containing the pathways/modules annotated to be linked to gene, obtained from 'load.gene2pathways' function in 'utils.R'
##' @param gene.id entrez gene id for the gene of interest
##' @param out.dir output directory for output files, data will be saved in two subdirectory in 'out.dir'
##' @return data files before and after applying Bonferroni correction and meta-analysis

## example:
# species <- 'human'
# data.files <- list.files('./data/output/human/merged/', full.names=T)
# sample.size <- load.sample.size(dir='./data/input/sample size/', species=species)
# r_mean <- load.r_mean(dir='./data/input/r_mean/', species=species)
# gene2pathways <- load.gene2pathways(dir='./data/utils data/gene2pathway/', species=species)
# gene.id <- 1
# out.dir <- './data/output/'

gmad_step2_gene <- function(data.files, sample.size, r_mean, species, gene2pathways, gene.id, out.dir){
  require(data.table)
  # dir.create(paste0(out.dir, 'GMAD_gene_preBonf/'), showWarnings=F)
  dir.create(paste0(out.dir, 'GMAD_gene/'), showWarnings=F)

  data.all <- NULL
  for(i in 1:length(data.files)){
    data.file <- data.files[i]
    data <- readRDS(data.file)
    class(data) <- 'data.frame'
    colnames(data)[1] <- 'path_id'
    tissue <- sapply(strsplit(data.file, split='//', fixed=TRUE), function(x) (x[2]))
    tissue <- gsub(pattern='_all_camera.RDS', replacement='', x=tissue, fixed=T)
    message('loading data from #', i, ' dataset: ', tissue)

    data.i <- data[, which(colnames(data) %in% c('path_id', gene.id))]
    if(is.null(ncol(data.i))) {
      message('dataset #', i, ': ', tissue, ' did not measure gene (ID: ', gene.id, ')')
      next
    }
    colnames(data.i)[2] <- tissue
    data.i <- data.table(data.i, key='path_id')
    if(is.null(data.all)){
      data.all <- data.i
    }else{
      data.all <- merge(data.all, data.i, all=TRUE)
    }
  }
  class(data.all) <- 'data.frame'

  # output
  # write.table(data.all, gzfile(paste0(out.dir, 'GMAD_gene_preBonf/', species,'--gene_',gene.id,'.gz')), quote=F, sep='\t', row.names=F, col.names=T)

  ### meta-analysis of -log10(p-values) from all datasets by taking weights of sample sizes and the average correlation coefficients for pathways in all datasets
  # apply Bonferroni correction on each datasets
  logp.thres <- -log10(0.05/nrow(data.all))
  data.sig <- data.matrix(data.all[, -which(colnames(data.all) %in% 'path_id')])
  rownames(data.sig) <- data.all$path_id
  data.sig[which((data.sig > -logp.thres) & (data.sig < logp.thres))] <- 0
  data.sig[which(data.sig >= logp.thres)] <- 1
  data.sig[which(data.sig <= -logp.thres)] <- -1
  data.sig <- data.frame(data.sig)

  # only keep datasets used
  # match the columns (datasets)
  sample.size$tissue <- make.names(sample.size$tissue)
  colnames(data.sig) <- make.names(colnames(data.sig))
  colnames(r_mean) <- make.names(colnames(r_mean))
  datasets <- intersect(intersect(sample.size$tissue, colnames(data.sig)), colnames(r_mean)) # get the datasets with all data

  sample.size.use <- sample.size[which(sample.size$tissue %in% datasets), ]
  r_mean.use <- r_mean[, which(colnames(r_mean) %in% datasets)]
  data.sig <- data.sig[, which(colnames(data.sig) %in% datasets)]
  sample.size.use <- sample.size.use[match(colnames(data.sig), sample.size.use$tissue), ]
  r_mean.use <- r_mean.use[, match(colnames(data.sig), colnames(r_mean.use))]

  # match the rows (pathways)
  r_mean.use <- r_mean.use[match(rownames(data.sig), rownames(r_mean.use)), ]

  if(nrow(sample.size.use) <= 1){
    print(paste0('gene: ', gene.id, 'is not measured in enough datasets'))
    next
  }

  # merge the gene correlations and sample size data
  r_mean_sample.size <- t(apply(r_mean.use, 1, function(x) x*sqrt(sample.size.use$size)))

  # compute the weighted.mean
  data.sig$mean.w <- rowSums(data.sig * r_mean_sample.size, na.rm=T) / rowSums(r_mean_sample.size, na.rm=T)
  data.sig$mean.w <- round(data.sig$mean.w, digits=3)

  # add whether the target gene is annotated to be in the pathway or not
  gene.pathways <- as.character(unlist(gene2pathways[gene.id]))
  data.all$pathway <- 0
  data.all$pathway[which(data.all$path_id %in% pathway.genes)] <- 1

  # output
  write.table(data.sig, gzfile(paste0(out.dir, 'GMAD_gene/', species,'--gene_',gene.id,'.gz')), quote=F, sep='\t', row.names=F, col.names=T)
}





################################################################################################################################################
######  analysis for one module/pathway vs all genes
##' @param data.files data files in RDS format for all datasets resulted from 'G-MAD_analysis_step1.R'
##' @param sample.size data.frame containing the sample size information for all datasets, obtained from 'load.sample.size' function in 'utils.R'
##' @param r_mean data.frame containing the average correlation coefficients for pathways in all datasets, obtained from 'load.r_mean' function in 'utils.R'
##' @param pathways list containing the pathway/module information, obtained from 'load.pathways' function in 'utils.R'
##' @param pathway.id id for the pathway/module of interest
##' @param out.dir output directory for output files, data will be saved in two subdirectory in 'out.dir'
##' @return data files before and after applying Bonferroni correction and meta-analysis

## example:
# species <- 'human'
# data.files <- list.files('~/Dropbox (Auwerx)/Hao Li/Gene set analysis/pgsa/data/output/human/merged/', full.names=T)
# data.files <- data.files[1:5]
# sample.size <- load.sample.size(dir='./data/input/sample size/', species=species)
# r_mean <- load.r_mean(dir='./data/input/r_mean/', species=species)
# pathways <- load.pathways(dir='./data/utils data/pathway/', species=species)
# pathway.id <- 'Reactome_R-HSA-191273'
# out.dir <- './data/output/'
# gmad_step2_module(data.files, sample.size, r_mean, species, pathway.id, out.dir)

gmad_step2_module <- function(data.files, sample.size, r_mean, species, pathways, pathway.id, out.dir){
  require(data.table)
  dir.create(paste0(out.dir, 'GMAD_module_preBonf/'), showWarnings=F)
  dir.create(paste0(out.dir, 'GMAD_module/'), showWarnings=F)

  data.all <- NULL
  for(i in 1:length(data.files)){
    data.file <- data.files[i]
    data <- readRDS(data.file)
    class(data) <- 'data.frame'
    rownames(data) <- data[,1]
    data <- data[, -1]
    # colnames(data)[1] <- 'path_id'
    tissue <- sapply(strsplit(data.file, split='//', fixed=TRUE), function(x) (x[2]))
    tissue <- gsub(pattern='_all_camera.RDS', replacement='', x=tissue, fixed=T)
    message('loading data from #', i, ' dataset: ', tissue)

    data.i <- data.frame(t(data[grep(pathway.id, rownames(data)), ]))
    data.i$gene_id <- rownames(data.i)

    if(is.null(ncol(data.i))) {
      message('dataset #', i, ': ', tissue, ' did not measure pathway (ID: ', pathway.id, ')')
      next
    }
    colnames(data.i)[1] <- tissue
    data.i <- data.table(data.i, key='gene_id')
    if(is.null(data.all)){
      data.all <- data.i
    }else{
      data.all <- merge(data.all, data.i, all=TRUE)
    }
  }
  class(data.all) <- 'data.frame'

  # save data.all for M-MAD
  write.table(data.all, gzfile(paste0(out.dir, 'GMAD_module_preBonf/',species,'--module_',pathway.id,'.gz')), quote=F, sep='\t', row.names=F, col.names=T)


  ### meta-analysis of -log10(p-values) from all datasets by taking weights of sample sizes and the average correlation coefficients for pathways in all datasets
  # apply Bonferroni correction on each datasets
  logp.thres <- -log10(0.05/nrow(data.all))
  data.sig <- data.matrix(data.all[, -which(colnames(data.all) %in% 'gene_id')])
  rownames(data.sig) <- data.all$gene_id
  data.sig[which((data.sig > -logp.thres) & (data.sig < logp.thres))] <- 0
  data.sig[which(data.sig >= logp.thres)] <- 1
  data.sig[which(data.sig <= -logp.thres)] <- -1
  data.sig <- data.frame(data.sig)

  # only keep datasets used
  # match the columns (datasets)
  sample.size$tissue <- make.names(sample.size$tissue)
  colnames(data.sig) <- make.names(colnames(data.sig))
  colnames(r_mean) <- make.names(colnames(r_mean))
  datasets <- intersect(intersect(sample.size$tissue, colnames(data.sig)), colnames(r_mean)) # get the datasets with all data

  # sort by the sample size of data.files
  sample.size.use <- sample.size[which(sample.size$tissue %in% datasets),] # should only use the common ones in match

  # get the gene correlations for the pathway and matched data.files
  r_mean.select <- r_mean[which(rownames(r_mean) == pathway.id), -1]
  r_mean.use <- data.frame(tissue=colnames(r_mean.select), r.mean=as.numeric(r_mean.select[1,]), stringsAsFactors=F)

  sample.size.use <- plyr::join(sample.size.use, r_mean.use, by='tissue', type='inner')
  sample.size.use <- sample.size.use[complete.cases(sample.size.use),]

  if(nrow(sample.size.use) <= 1){
    print(paste0('Module: ', pathway.id, 'is not analyzed in enough datasets'))
    next
  }

  data.sig <- data.sig[,match(make.names(sample.size.use$tissue), make.names(colnames(data.sig)))]
  data.sig <- data.frame(data.sig)

  data.sig$mean.w <- apply(data.sig, 1, FUN=function(x) weighted.mean(x, w=sqrt(sample.size.use$size)*sample.size.use$r.mean, na.rm=T))
  data.sig$mean.w <- round(data.sig$mean.w, digits=3)

  # add whether the target pathway/module is known to be linked to the genes
  pathway.genes <- as.character(unlist(pathways[pathway.id]))
  data.sig$pathway <- 0
  data.sig$pathway[which(data.sig$gene_id %in% pathway.genes)] <- 1

  data.sig$gene_id <- rownames(data.sig)
  data.sig <- data.sig[, c(ncol(data.sig), 1:(ncol(data.sig)-1))]
  # output
  write.table(data.sig, gzfile(paste0(out.dir, 'GMAD_module/',species,'--module_',pathway.id,'.gz')), quote=F, sep='\t', row.names=F, col.names=T)
}
