# GeneBridge
Systems analysis toolkit containing Gene-Module Association Determination (G-MAD) and Module-Module Association Determination (M-MAD) to identify the associations between genes and biological process or modules, as well as associations between modules.

For analyzed results for a given gene or module/gene set, please visit https://systems-genetics.org/genebridge.

## Citation
Li, H. et al. Identifying gene function and module connections by the integration of multi-species expression compendia. In revision. Genome Research 2019. https://genome.cshlp.org/content/early/2019/11/21/gr.251983.119.abstract.


## Dependencies

#### Probabilistic Estimation of Expression Residuals (PEER).

Installation: https://github.com/PMBio/peer/wiki/Installation-instructions.

Reference: Stegle O. et al. Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. Nat Protoc. 2012. https://www.nature.com/articles/nprot.2011.457

#### Install other dependent packages `data.table`, `plyr`, `parallel`, `caret`, `limma`, `impute`
```R
    install.packages('data.table')
    install.packages('plyr')
    install.packages('parallel')
    install.packages('caret')

    if (!requireNamespace('BiocManager', quietly=TRUE))
        install.packages('BiocManager')
    BiocManager::install('limma', version='3.8')
    BiocManager::install('impute')
```    
Reference: `limma` Ritchie ME. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015. https://academic.oup.com/nar/article/43/7/e47/2414268.

## Usage
```R
source('utils.R')
```    

#### Step 1: remove known and hidden covariates and get the expression residuals using PEER
```R
source('peer.R')

# load gene id to symbol information
mat.aligner <- read.table('./data/utils data/gene info/human_gene_info.txt', header=T, sep='\t', quote='')
mat.aligner <- mat.aligner[, c('gene_id', 'symbol')]
colnames(mat.aligner) <- c('id','gene')

data.dir <- './data/input/expr/'
expr.files <- list.files(data.dir, pattern='_expr.gct', full.names=T)
cov.files <- list.files(data.dir, pattern='_cov.txt', full.names=T)

# use GSE77688 as the example
tissue <- 'GSE77688'
expr.file <- grep(tissue, expr.files, value=T)
cov.file <- grep(tissue, cov.files, value=T)

if(length(cov.file) == 1){
  covs <- read.table(cov.file, header=T, row.names=1, sep='\t', quote='')
}else{
  covs <- NULL
}

# load data
mat <- read.table(expr.file, header=T, row.names=1, sep='\t', quote='')
mat <- mat[, -1]

# apply PEER
peer.exp <- peer.residual(mat, mat.aligner, covs, n.iteration=1000)

# output
out.file <- gsub(pattern='_expr.gct', replacement='_expr.peer.gct', x=expr.file)
out.file <- gsub(pattern='/expr/', replacement='/peer/', x=out.file)
write.table(peer.exp, out.file, quote=F, sep='\t', row.names=F, col.names=T)
```    

### load utils data
```R
species <- 'human'

# pathways: list containing the pathway/module information, obtained from 'load.pathways' function in 'utils.R'
pathways <- load.pathways(dir='./data/utils data/pathway/', species=species)
# pathway.pos: module position in the x axis (used for plotting)
pathway.pos <- load.pathway.pos(dir='./data/input/module position/', species=species)
# pathway.names: module names
pathway.names <- load.pathway.names(dir='./data/utils data/pathway/', species=species)

# gene2pathway: list containing the pathways/modules annotated to be linked to gene, obtained from 'load.gene2pathways' function in 'utils.R'
gene2pathways <- load.gene2pathways(dir='./data/utils data/gene2pathway/', species=species)

# gene.pos: gene position in the chromosomes (used for plotting)
gene.pos <- load.gene.pos(dir='./data/utils data/gene position/', species=species)

# sample.size: data.frame containing the sample size information for all datasets, obtained from 'load.sample.size' function in 'utils.R'
sample.size <- load.sample.size(dir='./data/input/sample size/', species=species)

# r_mean: data.frame containing the average correlation coefficients for pathways in all datasets, obtained from 'load.r_mean' function in 'utils.R'
r_mean <- load.r_mean(dir='./data/input/r_mean/', species=species)

# output directory
out.dir <- './data/output/'
dir.create(out.dir)

```    



#### Step 2: First step of G-MAD
```R
source('G-MAD_analysis_step1.R')

## load data
# data.files: data files after PEER correction
data.files <- list.files('./data/input/peer/', full.names=T)
# create output directory for output files
dir.create(paste0(out.dir,'GMAD_step1/'))

for(data.file in data.files){
  print(data.file)
  tissue <- sapply(strsplit(data.file, split='//', fixed=TRUE), function(x) (x[2]))
  tissue <- gsub(pattern='_expr.peer.gct', replacement='', x=tissue, fixed=T)

  results <- gmad_step1(data.file, pathways, num.cores=4) # change num.cores to enable parallel computing
  saveRDS(results, paste0(out.dir,'GMAD_step1/', tissue, '.RDS'))
}
```    


#### Step 3: Second step of G-MAD
```R
source('G-MAD_analysis_step2.R')

### G-MAD analysis for one gene against all modules

# data.files: data files in RDS format for all datasets resulted from G-MAD_step1 function in 'G-MAD_analysis_step1.R'
data.files <- list.files('./data/output/GMAD_step1/', pattern='RDS', full.names=T)      
# create output directory for output files, data will be saved in the subdirectory ('GMAD_gene') of out.dir
dir.create(paste0(out.dir,'GMAD_gene/'))

# gene.id: entrez gene id for the gene of interest
gene.id <- 1                

# run the gmad_step2_gene function
gmad_step2_gene(data.files, sample.size, r_mean, species, gene2pathways, gene.id, out.dir)


## make Manhattan plots for one gene against all modules
source('G-MAD_plot.R')

# load data from gmad_step2_gene function
data <- read.table('./data/output/GMAD_gene/human--gene_1652.gz', header=T, sep='\t', quote='')

# make plot
gmad_gene_plot(data, pathway.pos, pathway.names, species, thres=0.268, show.name=T)
```    


```R
### G-MAD analysis for one module against all genes
source('G-MAD_analysis_step2.R')

# data.files: data files in RDS format for all datasets resulted from G-MAD_step1 function in 'G-MAD_analysis_step1.R'
data.files <- list.files('./data/output/GMAD_step1/', pattern='RDS', full.names=T)
# create output directory for output files, data will be saved in two subdirectories ('GMAD_module_preBonf' and 'GMAD_module') of out.dir
dir.create(paste0(out.dir,'GMAD_module_preBonf/'))
dir.create(paste0(out.dir,'GMAD_module/'))

# pathway.id: id for the pathway/module of interest
pathway.id <- 'Reactome_R-HSA-611105'  

## run the gmad_step2_module function
gmad_step2_module(data.files, sample.size, r_mean, species, pathways, pathway.id, out.dir)


## make Manhattan plots for one module against all genes
source('G-MAD_plot.R')

# load data from gmad_step2_module function
data <- read.table('./data/output/GMAD_module/human--module_Reactome_R-HSA-191273.gz', header=T, sep='\t', quote='')

# make plot
gmad_module_plot(data, gene.pos, species, thres=0.268, show.name=T)
```    



#### Step 4. M-MAD
```R
source('M-MAD_analysis.R')

# data.files: gz files in 'GMAD_module_preBonf' subfolder for all modules resulted from gmad_step2_module function in 'G-MAD_analysis_step2.R'
data.files <- list.files('./data/output/GMAD_module_preBonf/', pattern='gz', full.names=T)
# create output directory for output files, data will be saved in the subdirectory ('MMAD_module') of out.dir
dir.create(paste0(out.dir, 'MMAD_module/'))

# pathway.id: id for the pathway/module of interest (use one dataset as the example )
pathway.id <- 'Reactome_R-HSA-611105'
data.file <- grep(pathway.id, data.files, value=T)
print(data.file)

# load data
data <- read.table(data.file, header=T, sep='\t', quote='', row.names=1)

# apply mmad_step1 and mmad_step2 functions
result <- mmad_step1(data, pathways)
result.sig <- mmad_step2(result, pathways, sample.size, r_mean, pathway.id)

# export
write.table(result.sig, gzfile(paste0(out.dir,'MMAD_module/',species,'--module_',pathway.id,'.gz')), quote=F, sep='\t', row.names=F, col.names=T)


## make Manhattan plots for one module against all modules
source('M-MAD_plot.R')

# load data
data <- read.table('./data/output/MMAD_module/human--module_GOBP_GO:0008610.gz', header=T, sep='\t', quote='')

# make plot
mmad_module_plot(data, pathway.pos, pathway.names, species, thres=0.268, show.name=T)
```

