#' dimarDoImputations
#'
#' @description Performs imputations
#' @importFrom foreach %dopar%
#' @return Imputed data set(s)
#' @param mtx Quantitative matrix
#' @param method Imputation method(s)
#' @param lib R packages of imputation methods (to be loaded in parallel loop 'foreach')
#' @export dimarDoImputations
#' @examples
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' Imputations <- dimarDoImputations(mtx, c('impSeqRob', 'ppca', 'imputePCA'))

dimarDoImputations <- function(mtx, method=NULL, lib=NULL) {
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
      lib[m] <- dimarGetLib(method[m])
    }
  }

  if (length(dim(mtx))>2) {
    # load parallel package
    doParallel::registerDoParallel()
  }
  Imp <- list()
  # loop over imputation methods
  for (m in 1:length(method)) {
    if (length(dim(mtx))>2) {
      # parallelize over simulated patterns (3rd dimension of mtx)
      I <- foreach::foreach(i=1:dim(mtx)[3], .combine='cbind', .packages=lib[m], .export=c("dimarDoImputationsR")) %dopar% {
        dimarDoImputationsR(mtx[,,i],method[m], lib[m])
      }
      if (!is.null(I)) {
        Imp[[method[m]]] <- array(as.matrix(I),dim=c(dim(mtx)[1],dim(mtx)[2],dim(mtx)[3]))
      } else { warning('DoImputations.R: Imputation of method ',method[m],' not feasible.')}
    } else {
      I <- dimarDoImputationsR(mtx,method[m],lib[m])
      if (!is.null(I)) {
        Imp[[method[m]]] <- as.matrix(I[,1:dim(mtx)[2]])
      } else { warning('DoImputations.R: Imputation of method ',method[m],' not feasible.')}
    }
    print(paste('Imputation with',method[m],'done.'))
  }
  # delete methods which have NA in imputation
  if (any(is.na(Imp))){
    naid <- which(is.na(Imp), arr.ind=T)
    Imp <- Imp[[!unique(naid[,4])]] # 4th dim are algorithms
  }
  return(Imp)
}

