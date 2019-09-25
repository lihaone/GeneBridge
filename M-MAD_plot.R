#!/usr/bin/env Rscript
# script to make MMAD plot for one pathway/module

##' @param data MMAD results for one pathway/module      # data.frame of MMAD results for the given pathway against all pathways (arranged in rows), should contain columns of pathway id ('path_id') and final MMAS ('mean.w')
##' @param pathway.pos pathway position information, obtained from 'load.pathway.pos' function
##' @param pathway.names pathway names information, obtained from 'load.pathway.names' function
##' @param species species for analysis
##' @param thres threshold for plotting. default as 0.268.
##' @return MMAD plot for one pathway

# # example:
# species <- 'human'
# data <- read.table('./data/output/MMAD_module/human--GOBP_GO:0008610.gz', header=T, sep='\t', quote='')
# pathway.pos <- load.pathway.pos(dir='./data/input/module position/', species=species)
# pathway.names <- load.pathway.names(dir='./data/utils data/pathway/', species=species)
#
# mmad_module_plot(data, pathway.pos, pathway.names, species, thres=0.268, show.name=T)

mmad_module_plot <- function(data, pathway.pos, pathway.names, species, thres=0.268, show.name=T){
  data <- data[, c('path_id','mean.w')]
  data <- data[order(data$mean.w, decreasing=T), ]
  data$color <- 'grey'
  data[which(abs(data$mean.w) >= thres), 'color'] <- 'black'
  data$size <- thres
  data[which(abs(data$mean.w) >= 0.1), 'size'] <- 3*abs(data[which(abs(data$mean.w) >= 0.1), 'mean.w'])

  library(plyr)
  # obtain the positions of pathways for plotting
  data <- join(data, pathway.pos, by='path_id', type='inner')
  # obtain pathway names
  data <- join(data, pathway.names, by='path_id', type='inner')

  # create plot
  ymax <- max(max(data$mean.w, na.rm=T), 0.4)
  ymin <- min(min(data$mean.w, na.rm=T), -0.4)
  plot(data$pos, data$mean.w, pch=19, cex=data$size, ylim=c(ymin,ymax), col=data$color,
       xlab='Modules', ylab='MMAS', xaxt='n', yaxt='n',
       cex.lab=1.5, cex.axis=1.2)
  axis(side=2, las=2)

  abline(h=0, col='black', lwd=1)
  abline(h=c(thres, -thres), col='red', lty=2, lwd=2)
  if(show.name & (length(which(abs(data$mean.w) >= thres)) != 0)){
    text(x=data[which(abs(data$mean.w) >= thres),'pos'],
         y=data[which(abs(data$mean.w) >= thres),'mean.w'],
         labels=data[which(abs(data$mean.w) >= thres),'path_name'])
  }
}
