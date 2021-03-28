#' dimarDoOptimalImputation
#'
#' @description Performs the optimal imputation to quantitative matrix.
#' @return Quantitative matrix on which optimal imputation method was applied.
#' @param mtx Quantitative matrix
#' @param method Imputation method
#' @param lib Library
#' @export dimarDoOptimalImputation
#' @examples Sample example to demonstrate the function
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

