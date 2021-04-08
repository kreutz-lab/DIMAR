# DIMAR
Data-driven selection of an imputation algorithm in R.

# Installation

For installation of DIMAR run:
```
install.packages("https://cran.r-project.org/src/contrib/Archive/imputation/imputation_1.3.tar.gz", repos=NULL, type='source')

devtools::install_github("cran/DMwR")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(pkgs = c("pcaMethods", "impute", "SummarizedExperiment"))

install.packages("devtools")
devtools::install_github("kreutz-lab/DIMAR")
```

# Examples:

MaxQuant example: DIMA applied to Reimann et al. (2020):
```
library(DIMAR)
filename <- "proteinGroups_PXD008893.txt"
filepath <- system.file("extdata", filename, package = "DIMAR")
Imp <- DIMAR::dimar(filepath,pattern='Intensity',group=c('PKB','PKC')
```
