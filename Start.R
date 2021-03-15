Sys.setenv(LANG="en")
library("SummarizedExperiment")
library("stats")
source("readInMaxQuant.R")
source("GetDesign.R")

filename <- "proteinGroups.txt"
out <- maxQuantToMatrix(filename)
out <- out[!rowData(out)[["ixs"]],] 

# extract data and feature annotation
mtx <- assays(out)[["data"]]
mtx <- mtx[1:10,1:5]
colnames(mtx) <- gsub("LFQ intensity","",colnames(mtx))
namtx <- is.na(mtx)

# Subsample indices
nsub <- ceiling(dim(mtx)[1]/10)
npersub <- ceiling(dim(mtx)[1]/nsub)
indrand <- sample(1:dim(mtx)[1],dim(mtx)[1])

for (i in 1:nsub) {
  ind <- indrand[npersub*(i-1)+1:npersub*i]
  X <- GetDesign(mtx[ind,])
  y <- as.numeric(namtx[ind,])
  
  #fit <- glm.fit(X,y,family=binomial(),weights=rep(1,dim(X)[1]))
  fit <- glm.fit(X[[1]],y, family=binomial())
  coef <- coefficients(fit)
}

coefmean <- c(mean(coef[X[[2]]~=3],na.rm=TRUE),coef[X[[2]]==3]) # row coefficients stay same, intensity/column coefficients are set to mean over nsub (for loop)
