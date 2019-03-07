# GeneBridge
Systems analysis toolkit containing Gene-Module Association Determination (G-MAD) and Module-Module Association Determination (M-MAD) to identify the associations between genes and biological process or modules, as well as associations between modules. 

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
