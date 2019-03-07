# !/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)

# script to analyze pathway interactions using data in paste0(species.i,"/merged/") folder
library(data.table)
library(limma)
options(stringsAsFactors = F)
set.seed(666)

main.dir <- "./"
setwd(main.dir)
species.i <- "human"

# get coding gene list
gene.info.dir <- 'data/gene/'
gene.info.files <- list.files(gene.info.dir, pattern = 'txt', full.names = T)
gene.info.file <- grep(species.i, gene.info.files, value = T)
gb <- read.table(gene.info.file, sep = '\t', header = T, quote = '', stringsAsFactors = F)
# remove ensemble ID for now
gb <- gb[which(gb$type == 'protein-coding'), -1]

# pathway size information
pathway.dir <- "data/pathway/"
pathway.files <- list.files(pathway.dir, pattern = "data", full.names = TRUE)
pathway.file <- grep(pattern = species.i, x = pathway.files, value = TRUE)
pathways <- readRDS(pathway.file)
# filter out gene sets with gene numbers less than 15 and larger than 1000
pathways <- pathways[which(lengths(pathways) >= 15)]
pathways <- pathways[which(lengths(pathways) <= 1000)]

path.name.files <- list.files(pathway.dir, pattern = "name", full.names = TRUE)
path.name.file <- grep(pattern = species.i, x = path.name.files, value = TRUE)
pathway.names <- readRDS(path.name.file)


data.dir <- paste0("data/example/output/",species.i,"/merged/")
data.files <- list.files(data.dir, pattern = ".RDS", full.names = T)

out.dir <- gsub(pattern = '/merged/', replacement = '/interpathway/', data.dir)
dir.create(out.dir)
dir.create(paste0(out.dir, "merged/"))
dir.create(paste0(out.dir, "separated/"))
dir.create(paste0(out.dir, "pathway/"))
dir.create(paste0(out.dir, "pathway_sign/"))

sample.size <- read.table(paste0("data/sample size/",species.i,"_sample_size.txt"), header=T, sep = "\t", quote = "", stringsAsFactors=F)
sample.size$tissue <- gsub(pattern = " ", replacement = "_", x = sample.size$tissue)
sample.size <- sample.size[order(sample.size$size, decreasing = F),]
# use data.files with sample.size larger than 80
sample.size <- sample.size[which(sample.size$size >= 80),]

# check for which datasets have been processed already !!!
out.files <- list.files(paste0(out.dir, "merged/"), pattern = 'RDS', full.names = F)
out.tissues <- gsub(pattern = '.RDS', replacement = '', x = out.files, fixed = T)

if((length(out.tissues) != 0) & (length(grep(paste0(out.tissues, collapse = '|'), data.files)) != 0)){
  data.files <- data.files[-grep(paste0(out.tissues, collapse = '|'), data.files)]
}

print(data.files)



for(data.file.i in 1:length(data.files)){
  data.file <- data.files[data.file.i]
  print(c(data.file.i, data.file))
  # data.file <- grep("GSE31705", data.files, value = T)

  dataset.i <- strsplit(data.file, split = "//", fixed = T)[[1]][2]
  dataset.i <- gsub(pattern = "_all_camera.RDS", replacement = "", x = dataset.i)
  data <- readRDS(data.file)
  class(data) <- 'data.frame'

  result.all <- mmad_step1(data, pathways)
  
  # save result.all
  saveRDS(result.all, paste0(out.dir, "merged/", dataset.i, ".RDS"))
}
