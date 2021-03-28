# DIMAR
Data-driven selction of an imputation algorithm in R.

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
mtx <- maxQuantToMatrix(filepath)

coef <- LearnPattern(mtx)
ref <- ConstructReferenceData(mtx)
sim <- AssignPattern(ref,coef,mtx)

Imputations <- DoImputations(sim,c('impSeqRob','ppca','imputePCA'))
Performance <- EvaluatePerformance(Imputations,ref,sim,'RMSE',TRUE)
Imp <- DoOptimalImputation(mtx,rownames(Performance))
write.csv(Imp, file=paste0("Imp_",filename))
```
