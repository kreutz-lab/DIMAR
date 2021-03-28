#' dimarDoImputations
#'
#' @description Performs imputations
#' @importFrom foreach %dopar%
#' @return Imputated dataset
#' @param mtx Quantitative matrix
#' @param method Imputation method
#' @param lib Library
#' @export dimarDoImputations
#' @examples Sample example to demonstrate the function
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

    # Initialize imputation array
    Imp <- array(NA,c(dim(mtx),length(method)))
  } else {
    Imp <- array(NA,c(dim(mtx), 1, length(method)))
  }

  # loop over imputation methods
  for (m in 1:length(method)) {
    if (length(dim(mtx))>2) {
      # parallelize over simulated patterns (3rd dimension of mtx)
      I <- foreach(i=1:dim(mtx)[3], .combine='cbind', .packages=lib[m], .export=c("dimarDoImputationsR")) %dopar% {
        dimarDoImputationsR(mtx[,,i],method[m], lib[m])
      }
      if (!is.null(I)) {
        Imp[,,,m] <- as.matrix(I[, 1:dim(mtx)[2]])
      } else { warning('DoImputations.R: Imputation of method ',method[m],' not feasible.')}
    } else {
      I <- dimarDoImputationsR(mtx,method[m],lib[m])
      if (!is.null(I)) {
        Imp[, , 1, m] <- as.matrix(I[,1:dim(mtx)[2]])
      } else { warning('DoImputations.R: Imputation of method ',method[m],' not feasible.')}
    }
    print(paste('Imputation with',method[m],'done.'))
  }
  # delete methods which have NA in imputation
  if (any(is.na(Imp))){
    naid <- which(is.na(Imp), arr.ind=T)
    Imp <- Imp[, , , !unique(id[,4])] # 4th dim are algorithms
    method <- method[!unique(id[,4])]
  }

  Imp <- list(Imp,method)
  names(Imp) <- c('Imp','method')
  return(Imp)
}

