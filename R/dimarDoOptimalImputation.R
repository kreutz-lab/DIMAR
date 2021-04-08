#' dimarDoOptimalImputation
#'
#' @description Performs imputation of the quantitative matrix with the given method (optimally selected by DIMA)
#' @return Imputed matrix
#' @param mtx Quantitative matrix
#' @param method Optimal imputation method
#' @param lib R package of imputation method
#' @export dimarDoOptimalImputation
#' @examples 
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' Imp <- dimarDoOptimalImputation(mtx, 'impSeq')

dimarDoOptimalImputation <- function(mtx, method, lib=NULL) {
m <- 0
Imp <- NULL
while (is.null(Imp) || any(is.na(Imp))) {
  m <- m+1
  eval(parse(text=paste('require(', dimarGetLib(method[m]),')')))
  Imp <- dimarDoImputationsR(mtx, method[m], dimarGetLib(method[m]))
}
Imp <- as.matrix(Imp[, 1:dim(mtx)[2]])
print(paste('Imputation of input data with algorithm', method[m], 'is performed.'))
return(Imp)
}

