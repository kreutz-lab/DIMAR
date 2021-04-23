# DIMAR
Data-driven selection of an imputation algorithm in R. 

For further inforamtion refer to the publication [Egert et al. (2020)](https://www.biorxiv.org/content/10.1101/2020.10.13.323618v1) or the [DIMAR wiki](https://github.com/kreutz-lab/DIMAR/wiki).

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


# Examples

DIMA can take a numeric matrix or the file path to a MaxQuant ProteinGroups file as an input. The data is reduced to the columns which include `pattern` in their sample names. The imputation algorithms can be defined by the user, by default the nine most frequently selected algorithms of 142 Pride data sets are applied.

MaxQuant file path as input (Example taken from (Reimann et al. (2020))):
```
library(DIMAR)
filename <- "Test1.txt"
filepath <- system.file("extdata", filename, package = "DIMAR")
Imp <- DIMAR::dimar(mtx = filepath, pattern = 'Intensity', group = c('PKB','PKC'))
```

Matrix as input 
(Example taken from Khoonsari, Payam Emami, et al. "Analysis of the cerebrospinal fluid proteome in Alzheimer's disease." PloS one 11.3 (2016): e0150672.): 

```
library(DIMAR)
library(openxlsx)
filename2 <- "Test2.xlsx"
filepath2 <- system.file("extdata", filename2, package = "DIMAR")
df <- openxlsx::read.xlsx(filepath2, sheet="Sheet1", startRow = 2) 
row.names(df) <- paste(c(1:nrow(df)), df$`Protein.(Uniprot.ID)`, sep = "_") 
mtx <- as.matrix(df[, grepl("^AD\\d|^C\\d", names(df))])
Imp2 <- DIMAR::dimar(mtx = mtx, pattern = NULL)
```
Same example with defining the imputation algorithms:
```
Imp2 <- DIMAR::dimar(mtx = mtx, pattern = "^AD\\d|^C\\d", methods = c('impSeqRob','impSeq','missForest','imputePCA','ppca','bpca'))
```
