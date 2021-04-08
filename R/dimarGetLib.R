#' dimarGetLib
#'
#' @description Gets R library of the imputation method
#' @return R library of entered imputation method
#' @param method Imputation method string
#' @export dimarGetLib
#' @examples
#' lib <- dimarGetLib('impSeq')

dimarGetLib <- function(method) {
  if (method == 'midastouch' || method == 'rf' || method == 'mean' || method == 'norm' ||
      method == 'ri' || method == 'pmm' || method == 'sample' || method == 'cart') {
    lib <- 'mice'
  } else if (method == 'knn' || method == 'impute.knn') {
    lib <- 'impute'
  } else if (method == 'impnorm' || method == 'imp.norm') {
    lib <- 'norm'
  } else if (method == 'Amelia' || method == 'amelia') {
    lib <- 'Amelia'
  } else if (method == 'regression' || method == 'aregImpute') {
    lib <- 'Hmisc'
  } else if (method == 'ppca' || method == 'bpca' || method == 'nipals' ||
             method == 'nlpca' || method == 'svd' || method == 'svdImpute') {
    lib <- 'pcaMethods'
  } else if (method == 'MinDet' || method == 'KNN' || method == 'MinProb' || method == 'QRILC') {
    lib <- 'imputeLCMD'
  } else if (method == 'SVTApproxImpute' || method == 'SVTImpute' || method == 'SVDImpute' ||
             method == 'kNNImpute' || method == 'lmImpute') {
    lib <- 'imputation'
  } else if (method == 'missForest') {
    lib <- 'missForest'
  } else if (method == 'softImpute') {
    lib <- 'softImpute'
  } else if (method == 'irmi') {
    lib <- 'VIM'
  } else if (method == 'Seq' || method == 'SeqRob' || method == 'impSeq' || method == 'impSeqRob') {
    lib <- 'rrcovNA'
  } else if (method == 'MIPCA' || method == 'imputePCA') {
    lib <- 'missMDA'
  } else if (method == 'mi') {
    lib <- 'mi'
  } else if (method == 'knnImputation') {
    lib <- 'DMwR'
  } else if (method == 'GMSimpute' || method == 'GMSLasso') {
    lib <- 'GMSimpute'
  } if (!exists('lib')) {
    warning(paste('dimarGetLib.R:',method,'is not implemented in DIMA and ignored for imputation. Check your spelling, or add respective Rcode to DoImputationsR.R and its library to GetLib.R.'))
  }
return(lib)
}
