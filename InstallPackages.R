
packages <- c('R.matlab','amap','mice','norm','Amelia','Hmisc','imputeLCMD','missForest','softImpute','VIM','rrcovNA','missMDA','mi','DMwR','GMSimpute')

for (i in 1:length(packages)){
  install.packages(packages[i], dependencies=TRUE, repos='http://cran.rstudio.com/')
}

# pcaMethods
install.packages("BiocManager")
BiocManager::install("pcaMethods")
BiocManager::install("impute")

# impute & imputation
install.packages("https://cran.r-project.org/src/contrib/Archive/imputation/imputation_1.3.tar.gz", repos=NULL, type='source')
