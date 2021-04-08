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
DIMA can take a numeric matrix or the file path to a MaxQuant ProteinGroups file as an input.

Matrix as input (Example taken from Khoonsari, Payam Emami, et al. "Analysis of the cerebrospinal fluid proteome in Alzheimer's disease." PloS one 11.3 (2016): e0150672.): 

```
library(DIMAR)
library(openxlsx)
filename2 <- "pone.0150672.s013.xlsx"
filepath2 <- system.file("extdata", filename2, package = "DIMAR")
df <- openxlsx::read.xlsx(filepath2, sheet="OpenMS_spiked_in_normalized_dat", startRow = 2) 
row.names(df) <- paste(c(1:nrow(df)), df$`Protein.(Uniprot.ID)`, sep = "_") 
mtx <- as.matrix(df[, grepl("^AD\\d|^C\\d", names(df))])
# Dataset reduced for demonstration purposes
mtx <- mtx[1:400,]
Imp2 <- DIMAR::dimar(mtx = mtx, pattern = NULL)
```

MaxQuant file path as input (Example taken from (Reimann et al. (2020))):
```
library(DIMAR)
filename <- "proteinGroups_PXD008893.txt"
filepath <- system.file("extdata", filename, package = "DIMAR")
Imp <- DIMAR::dimar(mtx = filepath, pattern = 'Intensity', group = c('PKB','PKC'))
```


