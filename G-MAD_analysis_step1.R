#!/usr/bin/env Rscript
# first step of G-MAD
# run pathway analysis for all genes using camera in each dataset

##' @param arrayData data.frame of expression residual with genes as rows and individuals as columns, obtained from 'peer.R'
##' @param pathways list of pathways containing the gene entrez id for each pathway, obtained from 'load.pathways' function in 'utils.R'
##' @param num.cores number of cores for the analysis, default=1
##' @return results.all data.table of camera p values between genes and pathways with pathways as rows and genes as columns

# # example:
# species <- 'human'
# data.files <- list.files('./data/input/peer/', full.names=T)
# pathways <- load.pathways(dir='./data/utils data/pathway/', species=species)
# out.dir <- './data/output/GMAD_step1/'
# dir.create(out.dir)
#
# for(data.file in data.files){
#   print(data.file)
#   tissue <- sapply(strsplit(data.file, split='//', fixed=TRUE), function(x) (x[2]))
#   tissue <- gsub(pattern='_expr.peer.gct', replacement='', x=tissue, fixed=T)
#
#   results.all <- gmad_step1(data.file, pathways)
#   saveRDS(results.all, paste0(out.dir, tissue, '.RDS'))
# }

gmad_step1 <- function(data.file, pathways, num.cores=1){
  library(limma)
  library(data.table)
  set.seed(666)

  arrayData <- read.table(data.file, header=T, sep='\t', quote='', row.names=1)
  arrayData <- arrayData[, -1]
  arrayData <- as.matrix(arrayData)

  # generate the index vector, made using ids2indices
  index <- ids2indices(pathways, rownames(arrayData), remove.empty=F)
  all.gene.num <- nrow(arrayData)

  if(num.cores == 1){                     # use only one core
    message('1 core is being used for the calculation.')
    for(target.gene.i in 1:all.gene.num){
      target.gene <- rownames(arrayData)[target.gene.i]
      message(target.gene.i, '/',all.gene.num)
      gene.i <- as.numeric(arrayData[target.gene.i, ])

      if(sd(gene.i, na.rm=T) == 0) {
        message(paste(target.gene, 'contains unique values, proceed with next gene'))
        next
      }

      design <- model.matrix( ~ gene.i)
      results <- camera(arrayData, index, design, inter.gene.cor=0.01, allow.neg.cor=T, use.ranks=F, sort=F)
      results$PValue <- results$PValue / 2

      results$path_id <- rownames(results)
      results$logp <- round(-log10(results$PValue), digits=3)
      results[which(results$Direction == 'Down'),'logp'] <- -results[which(results$Direction == 'Down'),'logp']
      results.x <- results[,which(colnames(results) %in% c('path_id', 'logp'))]
      colnames(results.x)[which(colnames(results.x) %in% c('logp'))] <- target.gene
      results.x <- data.table(results.x, key='path_id')

      if(target.gene.i == 1){
        results.all <- results.x
      }else{
        results.all <- merge(results.all, results.x, all=TRUE)
      }
    }
  }else{                         # use multiple cores
    # create a temperary directory for intermediate files
    temp.dir <- 'temp.dir/'
    dir.create(temp.dir)

    library(parallel)
    if(num.cores > detectCores()){
      stop('this device only have ',detectCores(),' cores, please set the "num.cores" lower than ',detectCores())
    }
    cl <- makeCluster(num.cores, type='FORK', outfile='')
    message(num.cores, ' cores are being used for the calculation.')

    batch.from <- 1
    batch.to <- num.cores
    batch.length <- ceiling(all.gene.num / num.cores)

    parLapply(cl, seq(batch.from, batch.to, by=1), function(batch.i){
      tryCatch({
        batch.i=as.integer(batch.i)
        gene.nos <- ((batch.i-1)*batch.length+1): min((batch.i*batch.length),all.gene.num)

        for(i in 1:length(gene.nos)){
          message(paste0('part ',batch.i, ', ', i, '/',length(gene.nos)))
          target.gene.i <- gene.nos[i]
          target.gene <- rownames(arrayData)[target.gene.i]
          gene.i <- as.numeric(arrayData[target.gene.i, ])

          if(sd(gene.i, na.rm=T) == 0) {
            message(paste(target.gene, 'contains unique values, proceed with next gene'))
            next
          }

          design <- model.matrix( ~ gene.i)
          results <- camera(arrayData, index, design, inter.gene.cor=0.01, allow.neg.cor=T, use.ranks=F, sort=F)
          results$PValue <- results$PValue / 2      # convert two-sided p-values into one-sided pvalues

          results$path_id <- rownames(results)
          results$logp <- round(-log10(results$PValue), digits=3)
          results[which(results$Direction == 'Down'),'logp'] <- -results[which(results$Direction == 'Down'),'logp']
          results.x <- results[,which(colnames(results) %in% c('path_id', 'logp'))]
          colnames(results.x)[which(colnames(results.x) %in% c('logp'))] <- target.gene
          results.x <- data.table(results.x, key='path_id')

          if(i == 1){
            results.all <- results.x
          }else{
            results.all <- merge(results.all, results.x, all=TRUE)
          }
        }

        message('saving part ', batch.i)
        saveRDS(results.all, paste0(temp.dir, 'temp_', batch.i, '.RDS'))

      }, error=function(err){
        message(paste('ERROR: ', err))
      })
    })
    stopCluster(cl)

    # merge the intermediate files from temp.dir
    data.files <- list.files(temp.dir, pattern='RDS', full.names=T)
    message('merging all parts')
    results.all <- NULL
    for(data.file in data.files){
      results.x <- readRDS(data.file)
      if(is.null(results.all)){
        results.all <- results.x
      }else{
        results.all <- merge(results.all, results.x, all=TRUE)
      }
    }
    file.remove(data.files)
    file.remove(temp.dir)

  }
  return(results.all)
}
