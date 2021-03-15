# DIMAR
Data-driven selction of an imputation algorithm in R.

Janine Egert, Bettina Warscheid, and Clemens Kreutz. *DIMA: Data-driven selection of a suitable imputation algorithm.* bioRxiv, 2020. doi: 10.1101/2020.10.13.323618.

# Installation

DIMA performs and evaluates various imputation algorithms. Therefore, the respective R packages sould be properly installed:
```
packages <- c('Rtools','R.matlab','amap','mice','norm','Amelia','Hmisc','imputeLCMD','missForest','softImpute','VIM','rrcovNA','missMDA','mi','DMwR','GMSimpute')
for (i in 1:length(packages)){
   install.packages(packages[i], dependencies=TRUE, repos='http://cran.rstudio.com/')
}
install.packages("BiocManager")
BiocManager::install("pcaMethods")
BiocManager::install("impute")
install.packages("https://cran.r-project.org/src/contrib/Archive/imputation/imputation_1.3.tar.gz", repos=NULL, type='source')
```
