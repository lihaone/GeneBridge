# utils functions used in the analysis
source(stringsAsFactors = F)

# load pathway data
load.pathway <- function(dir = "./data/pathway/", species='human'){
  pathway.files <- list.files(dir, pattern="data", full.names=T)
  pathway.file <- grep(pattern=species, x=pathway.files, value=T)
  pathways <- readRDS(pathway.file)
  
  # filter out gene sets with gene numbers less than 15 and larger than 1000
  pathways <- pathways[which(lengths(pathways) >= 15)]
  pathways <- pathways[which(lengths(pathways) <= 1000)]
  
  return(pathways)
}

# load pathway names
load.pathway.names <- function(dir = "./data/pathway/", species='human'){
  path.name.files <- list.files(dir, pattern="name", full.names=T)
  path.name.file <- grep(pattern=species, x=path.name.files, value=T)
  pathway.names <- readRDS(path.name.file)
  
  return(pathway.names)
}

# load pathway position data
load.pathway.pos <- function(dir = "./data/position/", species='human'){
  pathway.pos.files <- list.files(dir, pattern='_pathway_jaccard_index.position.txt', full.names=T)
  pathway.pos.file <- grep(species, pathway.pos.files, value=T)
  pathway.pos <- read.table(pathway.pos.file, sep='\t', header=T, quote='')
  
  return(pathway.pos)
}



# load gene id data
load.gene.id <- function(dir = "./data/gene/", species='human'){
  gene.info.files <- list.files(gene.info.dir, pattern='txt', full.names=T)
  gene.info.file <- grep(species, gene.info.files, value=T)
  gene.id <- read.table(gene.info.file, sep='\t', header=T, quote='')
  gene.id <- gene.id[which(gene.id$type == 'protein-coding'), -1]
    
  return(gene.id)
}


# load gene position data
load.gene.pos <- function(dir = "./data/position/", species='human'){
  gene.pos.files <- list.files(dir, pattern='_gene_symbol_pos.txt', full.names=T)
  gene.pos.file <- grep(species, gene.pos.files, value=T)
  gene.pos <- read.table(gene.pos.file, sep='\t', header=T, quote='')

  return(gene.pos)
}


# load sample size data
load.sample.size <- function(dir = "./data/sample size/", species='human'){
  sample.size <- read.table(paste0(dir,species,"_sample_size.txt"), header=T, sep="\t", quote="")
  sample.size$tissue <- gsub(pattern=" ", replacement="_", x=sample.size$tissue)
  sample.size <- sample.size[order(sample.size$size, decreasing=F),]
  # keep datasets with over 80 samples
  sample.size <- sample.size[which(sample.size$size >= 80),]
  
  return(sample.size)
}


# load average correlation coefficients (r_mean)
load.r_mean <- function(dir = "./data/r_mean/", species='human'){
  r_mean <- readRDS(paste0(dir, species, "_gene_r_mean.RDS"))
  class(r_mean) <- 'data.frame'
  
  return(r_mean)
}

