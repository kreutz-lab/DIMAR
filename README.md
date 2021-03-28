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
filepath <- system.file("extdata", "TestData.txt", package = "DIMAR")
mtx <- dimarReadInMaxQuant(filepath)

coef <- dimarLearnPattern(mtx)
ref <- dimarConstructReferenceData(mtx)
sim <- dimarAssignPattern(ref,coef,mtx)

Imputations <- dimarDoImputations(sim,c('impSeqRob','ppca','imputePCA'))
Performance <- dimarEvaluatePerformance(Imputations,ref,sim,'RMSE',TRUE)
Imp <- dimarDoOptimalImputation(mtx,rownames(Performance))
write.csv(Imp, file=paste0("Imp_",filename))
```
