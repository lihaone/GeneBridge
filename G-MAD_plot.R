#!/usr/bin/env Rscript
# script to make GMAD plot for one pathway or one gene

##' @param data GMAD results for one pathway.      # data.frame of GMAD results for the given pathway against all genes (arranged in rows), should contain columns of gene id ('gene_id'), final GMAS ('mean.w'), and whether the gene-pathway connection is known ('pathway')
##' @param gene.pos gene position information, obtained from "load.gene.pos" function
##' @param species species for analysis
##' @param thres threshold for plotting. default as 0.268. 
##' @return GMAD plot for pathway

gmad_pathway_plot <- function(data, gene.pos, species, thres=0.268){
  # keep only genes with data available from at least 30% of datasets
  data <- data[which(rowMeans(!is.na(data[,2:(ncol(data)-2)])) > 0.3),]
  
  data1 <- data[,c('gene_id','pathway','mean.w')]
  data1 <- data1[order(data1$mean.w, decreasing = T),]
  
  # obtain the positions of pathways for plotting
  data1 <- plyr::join(data1, gene.pos[,c(1,2,3,6)], by = "gene_id")
  data1 <- data1[!is.na(data1$chr),]
  
  # set alternate color for nearby chromosomes
  data1$chr.num <- as.numeric(data1$chr)
  data1[which(toupper(data1$chr) %in% 'X'),'chr.num'] <- max(data1$chr.num, na.rm = T) + 1
  data1[which(toupper(data1$chr) %in% 'Y'),'chr.num'] <- max(data1$chr.num, na.rm = T) + 1
  data1[which(toupper(data1$chr) %in% c('M','MT')),'chr.num'] <- max(data1$chr.num, na.rm = T) + 1
  data1$color <- 'gray'
  data1[which((data1$chr.num %% 2) == 0),'color'] <- 'gray50'
  data1[which(abs(data1$mean.w) >= 0.3),'color'] <- "black"
  data1[which(data1$pathway ==1),'color'] <- 'red'
  data1$size <- 0.2
  data1[which(abs(data1$mean.w) >= 0.1),'size'] <- 2*abs(data1[which(abs(data1$mean.w) >= 0.1),'mean.w'])

  # get the x axis labels (chromosomes)
  Cat <- unique(data1$chr)
  pos.labels <- vector(length=length(Cat))
  names(pos.labels) <- Cat
  for (i in 1:length(Cat)){
    a <- Cat[i]
    b <- subset(data1, data1$chr==a)
    pos.labels[i] <- (max(b$pos) + min(b$pos)) / 2
    rm(a,b)
  }
  if(species %in% c('human','rat')){
    pos.labels <- pos.labels[which(names(pos.labels) %in% c(1,5,10,15,20,'X','Y','M'))]
  }else if(species == 'mouse'){
    pos.labels <- pos.labels[which(names(pos.labels) %in% c(1,5,10,15,19,'X','Y','M'))]
  }
  
  ymax <- max(max(data1$mean.w), 0.4)
  ymin <- min(min(data1$mean.w), -0.4)
  plot(data1$pos, data1$mean.w, pch=19, cex=data1$size, ylim=c(ymin,ymax), col=data1$color, xlab='Chromosomes', ylab='GMAS', xaxt='n', las=2,
       cex.lab=1.5, cex.axis=1.2, bty='n')
  box(lwd=2)
  abline(h=0, col='black', lwd=1)
  abline(h=c(thres, -thres), col='red', lty=2, lwd=2)
  if(length(which(abs(data1$mean.w) >= thres)) != 0){
    text(x=data1[which(abs(data1$mean.w) >= thres),'pos'], y=data1[which(abs(data1$mean.w) >= thres),'mean.w'], labels=data1[which(abs(data1$mean.w) >= thres),'symbol'])
  }
  text(pos.labels, par("usr")[3] - 0.05, srt = 0, adj = 1, labels = names(pos.labels), xpd = TRUE, cex=1.3)
  
}



##' @param data GMAD results for one gene      # data.frame of GMAD results for the given gene against all pathways (arranged in rows), should contain columns of pathway id ('path_id'), final GMAS ('mean.w'), and whether the gene-pathway connection is known ('pathway')
##' @param pathway.pos pathway position information, obtained from "load.pathway.pos" function
##' @param pathway.names pathway names information, obtained from "load.pathway.names" function
##' @param species species for analysis
##' @param thres threshold for plotting. default as 0.268. 
##' @return GMAD plot for one gene

gmad_gene_plot <- function(data, pathway.pos, pathway.names, species, thres=0.268){
  
  data1 <- data[,c('path_id','pathway','mean.w')]
  data1 <- data1[order(data1$mean.w, decreasing = T),]
  
  data1$color <- 'grey'
  data1[which(abs(data1$mean.w) >= thres),'color'] <- "black"
  data1[which(data1$pathway == 1),'color'] <- 'red'
  
  data1$size <- thres
  data1[which(abs(data1$mean.w) >= 0.1),'size'] <- 3*abs(data1[which(abs(data1$mean.w) >= 0.1),'mean.w'])
  
  library(plyr)
  # obtain the positions of pathways for plotting
  data1 <- join(data1, pathway.pos, by = "path_id")
  # obtain pathway names
  data1 <- join(data1, pathway.names, by = "path_id")
  
  ymax <- max(max(data1$mean.w), 0.4)
  ymin <- min(min(data1$mean.w), -0.4)
  plot(data1$pos, data1$mean.w, pch=19, cex=data1$size, ylim=c(ymin,ymax), col=data1$color, xlab='Modules', ylab='GMAS', xaxt='n', las=2, 
       cex.lab=1.5, cex.axis=1.2, bty='n')
  box(lwd=2)
  
  abline(h=0, col='black', lwd=1)
  abline(h=c(thres, -thres), col='red', lty=2, lwd=2)
  if(length(which(abs(data1$mean.w) >= thres)) != 0){
    text(x=data1[which(abs(data1$mean.w) >= thres),'pos'], y=data1[which(abs(data1$mean.w) >= thres),'mean.w'], labels=data1[which(abs(data1$mean.w) >= thres),'path_name'])
  }
  
}


