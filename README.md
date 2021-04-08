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

Example usage of package:
```
library(DIMAR)
filename <- "TestData.txt"
filepath <- system.file("extdata", filename, package = "DIMAR")
Imp <- DIMAR::dimar(filepath)
```
