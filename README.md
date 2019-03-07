# GeneBridge
Systems analysis toolkit containing Gene-Module Association Determination (G-MAD) and Module-Module Association Determination (M-MAD) to identify the associations between genes and biological process or modules, as well as associations between modules. 

For analyzed results for a given gene or module/gene set, please visit https://systems-genetics.org/genebridge. 

## Dependencies

PEER. Installation: https://github.com/PMBio/peer/wiki/Installation-instructions

Install other dependent packages `data.table`, `plyr`, `limma`
```R
    install.packages('data.table')
    install.packages('plyr')
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("limma", version = "3.8")
```    

### Citation
Li, H. et al. GeneBridge, a toolset to impute gene function by integration of multi-species expression compendia. Submitted. 
