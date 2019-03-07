#!/usr/bin/env Rscript
# script to make MMAD plot for one pathway/module

##' @param data MMAD results for one pathway/module      # data.frame of MMAD results for the given pathway against all pathways (arranged in rows), should contain columns of pathway id ('path_id') and final MMAS ('mean.w')
##' @param pathway.pos pathway position information, obtained from "load.pathway.pos" function
##' @param pathway.names pathway names information, obtained from "load.pathway.names" function
##' @param species species for analysis
##' @param thres threshold for plotting. default as 0.268. 
##' @return MMAD plot for one pathway

mmad_pathway_plot <- function(data, pathway.pos, pathway.names, species, thres=0.268){
  data1 <- data[,c('path_id','mean.w')]
  data1 <- data1[order(data1$mean.w, decreasing = T),]
  data1$color <- 'grey'
  data1[which(abs(data1$mean.w) >= thres),'color'] <- "black"
  data1$size <- thres
  data1[which(abs(data1$mean.w) >= 0.1),'size'] <- 3*abs(data1[which(abs(data1$mean.w) >= 0.1),'mean.w'])
  
  library(plyr)
  # obtain the positions of pathways for plotting
  data1 <- join(data1, pathway.pos, by = "path_id")
  # obtain pathway names
  data1 <- join(data1, pathway.names, by = "path_id")
  # obtain pathway similarity
  path.sim.i <- data.frame(path_id=rownames(path.sim), path.sim=path.sim[,which(colnames(path.sim) == pathway.id)])
  data1 <- join(data1, path.sim.i, by = "path_id")
  
  ymax <- max(max(data1$mean.w), 0.4)
  ymin <- min(min(data1$mean.w), -0.4)
  plot(data1$pos, data1$mean.w, pch=19, cex=data1$size, ylim=c(ymin,ymax), col=data1$color, xlab='Modules', ylab='MMAS', xaxt='n',yaxt='n', 
       cex.lab=1.5, cex.axis=1.2)
  axis(side=2, at=seq(from=ymin, to=ymax, by=bar.steps[step.pick]), las=2)
  
  abline(h=0, col='black', lwd=1)
  abline(h=c(thres, -thres), col='red', lty=2, lwd=2)
  if(length(which(abs(data1$mean.w) >= thres)) != 0){
    text(x=data1[which(abs(data1$mean.w) >= thres),'pos'], y=data1[which(abs(data1$mean.w) >= thres),'mean.w'], labels=data1[which(abs(data1$mean.w) >= thres),'path_name'])
  }
}



