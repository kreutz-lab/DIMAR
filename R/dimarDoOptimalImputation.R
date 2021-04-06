#' dimarDoOptimalImputation
#'
#' @description Performs imputation of the quantitative matrix with the given method (optimally selected by DIMA)
#' @return Imputed matrix
#' @param mtx Quantitative matrix
#' @param method Optimal imputation method
#' @param lib R package of imputation method
#' @export dimarDoOptimalImputation
#' @examples 
#' Imputations <- dimarDoImputations(mtx)
#' Performance <- dimarEvaluatePerformance(Imputations, ref, mtx)
#' Imp <- dimarDoOptimalImputation(mtx, rownames(Performance))

dimarDoOptimalImputation <- function(mtx, method, lib=NULL) {
m <- 0
Imp <- NULL
while (is.null(Imp) || any(is.na(Imp))) {
  m <- m+1
  eval(parse(text=paste('require(', dimarGetLib(method[m]),')')))
  Imp <- dimarDoImputationsR(mtx, method[m], dimarGetLib(method[m]))
}
Imp <- as.matrix(Imp[, 1:dim(mtx)[2]])
print(paste('Imputation with', method[m], 'done.'))
return(Imp)
}

