source("GetLib.R")
source("DoImputationsR.R")

DoImputations <- function(mtx,method=NULL,lib=NULL) {
  
  # Assign imputation methods
  if (is.null(method)) {
    method <- c('impSeqRob','impSeq','missForest','imputePCA','ppca','MIPCA','bpca','SVDImpute','kNNImpute','regression','aregImpute','softImpute','MinDet','amelia','SVTImpute','irmi','knn','QRILC','nipals','MinProb','rf','sample','pmm','svdImpute','norm','cart','midastouch','mean','ri')
  } else if (method[1]=='fast') {
    method <- c('impSeqRob','impSeq','missForest','imputePCA','ppca','MIPCA','bpca','SVDImpute','kNNImpute')
  }
  
  # Assign R package library
  if (is.null(lib)) {
    lib <- character(length(method))
    for (m in 1:length(method)) {
      lib[m] <- GetLib(method[m])
    }
  }
  
  if (length(dim(mtx))>2) {
    # load parallel package
    tryCatch({ require(doParallel)},
      warning=function(c) {install.packages("doParallel")
        require(doParallel) })
    registerDoParallel()
    
    # Initialize imputation array
    Imp <- array(NA,c(dim(mtx),length(method)))
  } else { Imp <- array(NA,c(dim(mtx),1,length(method))) }
  
  # loop over imputation methods
  for (m in 1:length(method)) {
    print(m)
    tryCatch({ require(lib[m],character.only=T)}, # character.only needed when loading 'character'
      warning=function(c) {install.packages(lib[m])
        require(lib[m],character.only=T) })
    print(method[m])
    if (length(dim(mtx))>2) {
      # parallelize over simulated patterns (3rd dimension of mtx)
      I <- foreach(i=1:dim(mtx)[3],.combine='cbind',.packages=lib[m]) %dopar% {
        DoImputationsR(mtx[,,i],method[m],lib[m])
      }
      if (!is.null(I)) {
        Imp[,,,m] <- I
      } else { warning('DoImputations.R: Imputation of method ',method[m],' not possible.')}
    } else {
      I <- DoImputationsR(mtx,method[m],lib[m])
      if (!is.null(I)) {
        Imp[,,1,m] <- as.matrix(I[,1:dim(mtx)[2]])
      } else { warning('DoImputations.R: Imputation of method ',method[m],' not possible.')}
    }
  }
  return(Imp)
}

