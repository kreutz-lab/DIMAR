Sys.setenv(LANG="en")
library("SummarizedExperiment")
library("stats")
source("readInMaxQuant.R")
source("ConstructDesignMatrix.R")
source("ConstructRegularizationMatrix.R")
source("LearnPattern.R")

filename <- "proteinGroups.txt"
out <- maxQuantToMatrix(filename)
out <- out[!rowData(out)[["ixs"]],] 

# extract data and feature annotation
mtx <- assays(out)[["data"]]
mtx <- mtx[1:10,1:5]
colnames(mtx) <- gsub("LFQ intensity","",colnames(mtx))

coef <- LearnPattern(mtx)