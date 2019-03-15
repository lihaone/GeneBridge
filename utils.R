# utils functions used in the analysis
# loading data.frame using read.table/read.csv
options(stringsAsFactors = F)

# pdf setting
pdf.options(family = 'ArialMT', useDingbats = F)

# cn (corner) function to get the corner of data.frame
cn <- function(x, n=5) { x[1:min(n,nrow(x)), 1:min(n,ncol(x))] }

# load pathway data
load.pathways <- function(dir = "./data/utils data/pathway/", species='human'){
  pathway.files <- list.files(dir, pattern="data", full.names=T)
  pathway.file <- grep(pattern=species, x=pathway.files, value=T)
  pathways <- readRDS(pathway.file)
  
  # filter out gene sets with gene numbers less than 15 and larger than 1000
  pathways <- pathways[which(lengths(pathways) >= 15)]
  pathways <- pathways[which(lengths(pathways) <= 1000)]
  
  return(pathways)
}

# load pathway names
load.pathway.names <- function(dir = "./data/utils data/pathway/", species='human'){
  path.name.files <- list.files(dir, pattern="name", full.names=T)
  path.name.file <- grep(pattern=species, x=path.name.files, value=T)
  pathway.names <- readRDS(path.name.file)
  # 
  return(pathway.names)
}

# load pathway position data
load.pathway.pos <- function(dir = "./data/input/module position/", species='human'){
  pathway.pos.files <- list.files(dir, pattern='_pathway_jaccard_index.position.txt', full.names=T)
  pathway.pos.file <- grep(species, pathway.pos.files, value=T)
  pathway.pos <- read.table(pathway.pos.file, sep='\t', header=T, quote='')
  
  return(pathway.pos)
}

# load gene2pathway data
load.gene2pathways <- function(dir = "./data/utils data/gene2pathway/", species='human'){
  pathway.files <- list.files(dir, pattern="_gene2path.RDS", full.names=T)
  pathway.file <- grep(pattern=species, x=pathway.files, value=T)
  gene2pathways <- readRDS(pathway.file)
  
  return(gene2pathways)
}


# load gene id data
load.gene.id <- function(dir = "./data/utils data/gene info/", species='human'){
  gene.info.files <- list.files(gene.info.dir, pattern='txt', full.names=T)
  gene.info.file <- grep(species, gene.info.files, value=T)
  gene.id <- read.table(gene.info.file, sep='\t', header=T, quote='')
  gene.id <- gene.id[which(gene.id$type == 'protein-coding'), -1]
    
  return(gene.id)
}


# load gene position data
load.gene.pos <- function(dir = "./data/utils data/gene position/", species='human'){
  gene.pos.files <- list.files(dir, pattern='_gene_symbol_pos.txt', full.names=T)
  gene.pos.file <- grep(species, gene.pos.files, value=T)
  gene.pos <- read.table(gene.pos.file, sep='\t', header=T, quote='')

  return(gene.pos)
}


# load sample size data
load.sample.size <- function(dir = "./data/input/sample size/", species='human'){
  sample.size <- read.table(paste0(dir,species,"_sample_size.txt"), header=T, sep="\t", quote="")
  sample.size$tissue <- gsub(pattern=" ", replacement="_", x=sample.size$tissue)
  sample.size <- sample.size[order(sample.size$size, decreasing=F),]
  # keep datasets with over 80 samples
  sample.size <- sample.size[which(sample.size$size >= 80),]
  
  return(sample.size)
}


# load average correlation coefficients (r_mean)
load.r_mean <- function(dir = "./data/input/r_mean/", species='human'){
  r_mean <- readRDS(paste0(dir, species, "_gene_r_mean.RDS"))
  class(r_mean) <- 'data.frame'
  rownames(r_mean) <- r_mean$pathway
  r_mean <- r_mean[, -1]
  
  return(r_mean)
}

