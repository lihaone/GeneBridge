#!/usr/bin/env Rscript
# script to make GMAD plot for one module/pathway or one gene

##' @param data GMAD results for one module/pathway obtained from "gmad_step2_module" function in "G-MAD_step2.R".
##'             Data.frame of GMAD results for the given pathway against all genes (arranged in rows),
##'             should contain columns of gene id ('gene_id'), final GMAS ('mean.w'), and whether the gene-module/pathway connection is known ('pathway')
##' @param gene.pos gene position information, obtained from "load.gene.pos" function
##' @param species species for analysis
##' @param thres threshold for plotting. default as 0.268.
##' @return GMAD plot for one module/pathway.

# # example:
# data <- read.table('./data/output/GMAD_module/human--module_Reactome_R-HSA-191273.gz', header=T, sep='\t', quote="")
# gene.pos <- load.gene.pos(dir="./data/utils data/gene position/", species='human')
# species <- 'human'
# 
# gmad_module_plot(data, gene.pos, species, thres=0.268)

gmad_module_plot <- function(data, gene.pos, species, thres=0.268){
  # keep only genes with data available from at least 30% of datasets
  data <- data[which(rowMeans(!is.na(data[, 2:(ncol(data)-2)])) > 0.3), ]

  data <- data[,c('gene_id','pathway','mean.w')]
  data <- data[order(data$mean.w, decreasing=T),]

  # obtain the positions of genes for plotting
  data <- plyr::join(data, gene.pos[, c(1,2,3,6)], by="gene_id", type='inner')
  data <- data[!is.na(data$chr),]

  # set alternate color for nearby chromosomes
  data$chr.num <- as.numeric(data$chr)
  data[which(toupper(data$chr) %in% 'X'),'chr.num'] <- max(data$chr.num, na.rm=T) + 1
  data[which(toupper(data$chr) %in% 'Y'),'chr.num'] <- max(data$chr.num, na.rm=T) + 1
  data[which(toupper(data$chr) %in% c('M','MT')),'chr.num'] <- max(data$chr.num, na.rm=T) + 1
  data$color <- 'gray'
  data[which((data$chr.num %% 2) == 0),'color'] <- 'gray50'
  data[which(abs(data$mean.w) >= 0.3),'color'] <- "black"
  data[which(data$pathway ==1),'color'] <- 'red'
  data$size <- 0.2
  data[which(abs(data$mean.w) >= 0.1),'size'] <- 2*abs(data[which(abs(data$mean.w) >= 0.1),'mean.w'])

  # get the x axis labels (chromosomes)
  Cat <- unique(data$chr)
  pos.labels <- vector(length=length(Cat))
  names(pos.labels) <- Cat
  for (i in 1:length(Cat)){
    a <- Cat[i]
    b <- subset(data, data$chr==a)
    pos.labels[i] <- (max(b$pos) + min(b$pos)) / 2
    rm(a,b)
  }
  if(species %in% c('human','rat')){
    pos.labels <- pos.labels[which(names(pos.labels) %in% c(1,5,10,15,20,'X','Y','M'))]
  }else if(species == 'mouse'){
    pos.labels <- pos.labels[which(names(pos.labels) %in% c(1,5,10,15,19,'X','Y','M'))]
  }

  ymax <- max(max(data$mean.w, na.rm=T), 0.4)
  ymin <- min(min(data$mean.w, na.rm=T), -0.4)
  plot(data$pos, data$mean.w, pch=19, cex=data$size, ylim=c(ymin,ymax), col=data$color, xlab='Chromosomes', ylab='GMAS', xaxt='n', las=2,
       cex.lab=1.5, cex.axis=1.2, bty='n')
  box(lwd=2)
  abline(h=0, col='black', lwd=1)
  abline(h=c(thres, -thres), col='red', lty=2, lwd=2)
  if(length(which(abs(data$mean.w) >= thres)) != 0){
    text(x=data[which(abs(data$mean.w) >= thres),'pos'], y=data[which(abs(data$mean.w) >= thres),'mean.w'], labels=data[which(abs(data$mean.w) >= thres),'symbol'])
  }
  text(pos.labels, par("usr")[3] - 0.05, srt=0, adj=1, labels=names(pos.labels), xpd=TRUE, cex=1.3)
}



##' @param data GMAD results for one gene obtained from "gmad_step2_gene" function in "G-MAD_step2.R".
##'             data.frame of GMAD results for the given gene against all pathways/modules (arranged in rows),
##'             should contain columns of pathway id ('path_id'), final GMAS ('mean.w'), and whether the gene-pathway/module connection is known ('pathway')
##' @param pathway.pos pathway/module position information, obtained from "load.pathway.pos" function
##' @param pathway.names pathway/module names information, obtained from "load.pathway.names" function
##' @param species species for analysis
##' @param thres threshold for plotting. default as 0.268.
##' @return GMAD plot for one gene

# # example:
# data <- read.table('./data/output/GMAD_gene/human--gene_1652.gz', header=T, sep='\t', quote="")
# pathway.pos <- load.pathway.pos(dir="./data/input/module position/", species='human')
# pathway.names <- load.pathway.names(dir="./data/utils data/pathway/", species='human')
# species <- 'human'
# 
# gmad_gene_plot(data, pathway.pos, pathway.names, species, thres=0.268)


gmad_gene_plot <- function(data, pathway.pos, pathway.names, species, thres=0.268){
  # keep only genes with data available from at least 30% of datasets
  data <- data[which(rowMeans(!is.na(data[,2:(ncol(data)-2)])) > 0.3), ]

  data <- data[,c('path_id','pathway','mean.w')]
  data <- data[order(data$mean.w, decreasing=T),]

  data$color <- 'grey'
  data[which(abs(data$mean.w) >= thres),'color'] <- "black"
  data[which(data$pathway == 1),'color'] <- 'red'

  data$size <- thres
  data[which(abs(data$mean.w) >= 0.1),'size'] <- 3*abs(data[which(abs(data$mean.w) >= 0.1),'mean.w'])

  library(plyr)
  # obtain the positions of pathways/modules for plotting
  data <- join(data, pathway.pos, by="path_id", type='inner')
  # obtain pathway names
  data <- join(data, pathway.names, by="path_id", type='inner')

  ymax <- max(max(data$mean.w, na.rm=T), 0.4)
  ymin <- min(min(data$mean.w, na.rm=T), -0.4)
  plot(data$pos, data$mean.w, pch=19, cex=data$size, ylim=c(ymin,ymax), col=data$color, xlab='Modules', ylab='GMAS', xaxt='n', las=2,
       cex.lab=1.5, cex.axis=1.2, bty='n')
  box(lwd=2)

  abline(h=0, col='black', lwd=1)
  abline(h=c(thres, -thres), col='red', lty=2, lwd=2)
  if(length(which(abs(data$mean.w) >= thres)) != 0){
    text(x=data[which(abs(data$mean.w) >= thres),'pos'], y=data[which(abs(data$mean.w) >= thres),'mean.w'], labels=data[which(abs(data$mean.w) >= thres),'path_name'])
  }
}
