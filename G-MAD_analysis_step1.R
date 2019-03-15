#!/usr/bin/env Rscript
# first step of G-MAD
# run pathway analysis for all genes using camera in each dataset

##' @param arrayData data.frame of expression residual with genes as rows and individuals as columns, obtained from "peer.R"
##' @param pathways list of pathways containing the gene entrez id for each pathway, obtained from "load.pathways" function in "utils.R"
##' @return results.all data.table of camera p values between genes and pathways with pathways as rows and genes as columns

# # example:
# data.files <- list.files('./data/input/peer/', full.names=T)
# pathways <- load.pathways(dir="./data/utils data/pathway/", species='human')
# out.dir <- './data/output/GMAD_step1/'
#
# for(data.file in data.files){
#   print(data.file)
#   tissue <- sapply(strsplit(data.file, split='//', fixed=TRUE), function(x) (x[2]))
#   tissue <- gsub(pattern="_expr.peer.gct", replacement="", x=tissue, fixed=T)
#
#   results.all <- gmad_step1(data.file, pathways)
#   saveRDS(results.all, paste0(out.dir, tissue, '_all_camera.RDS'))
# }

gmad_step1 <- function(data.file, pathways){
  require(limma)
  require(data.table)
  set.seed(666)

  arrayData <- read.table(data.file, header=T, sep='\t', quote="", row.names=1)
  arrayData <- as.matrix(arrayData)

  # generate the index vector, made using ids2indices
  index <- ids2indices(pathways, rownames(arrayData), remove.empty=F)

  all.gene.num <- nrow(arrayData)
  for(target.gene.i in 1:all.gene.num){
    target.gene <- rownames(arrayData)[target.gene.i]
    message(paste0(target.gene.i, '/',all.gene.num))
    gene.i <- as.numeric(arrayData[target.gene.i, ])

    if(sd(gene.i) == 0) {
      message(paste(target.gene, "contains unique values, proceed with next gene"))
      next
    }

    design <- model.matrix( ~ gene.i)
    results <- camera(arrayData, index, design, inter.gene.cor=0.01, allow.neg.cor=T, use.ranks=F, sort=F)
    results$PValue <- results$PValue / 2    

    results$path_id <- rownames(results)
    results$logp <- round(-log10(results$PValue), digits=3)
    results[which(results$Direction == "Down"),"logp"] <- -results[which(results$Direction == "Down"),"logp"]
    results.x <- results[,which(colnames(results) %in% c("path_id", "logp"))]
    colnames(results.x)[which(colnames(results.x) %in% c("logp"))] <- target.gene
    results.x <- data.table(results.x, key="path_id")

    if(target.gene.i == 1){
      results.all <- results.x
    }else{
      results.all <- merge(results.all, results.x, all=TRUE)
    }
  }
  return(results.all)
}
