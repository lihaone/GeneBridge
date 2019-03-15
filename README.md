# GeneBridge
Systems analysis toolkit containing Gene-Module Association Determination (G-MAD) and Module-Module Association Determination (M-MAD) to identify the associations between genes and biological process or modules, as well as associations between modules.

For analyzed results for a given gene or module/gene set, please visit https://systems-genetics.org/genebridge.

## Dependencies

Probabilistic Estimation of Expression Residuals (PEER). 

Installation: https://github.com/PMBio/peer/wiki/Installation-instructions. 

Reference: Stegle O. et al. Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. Nat Protoc. 2012 

Install other dependent packages `data.table`, `plyr`, `limma`
```R
    install.packages('data.table')
    install.packages('plyr')

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("limma", version = "3.8")
```    
Reference: 

limma: Ritchie ME. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015. 


### Citation
Li, H. et al. Identifying gene function and module connections by the integration of multi-species expression compendia. Submitted. 
