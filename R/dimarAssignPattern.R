#' dimarAssignPattern
#'
#' @description Assigns patterns of MVs to reference data set with a logistic regression model
#' @return Reference dataset with applied pattern
#' @param ref Reference dataset
#' @param coef Logistic regression coefficients of missing value pattern
#' @param mtx Quantitative matrix
#' @param npat Number of patterns
#' @export dimarAssignPattern
#' @examples
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' coef <- dimarLearnPattern(mtx)
#' ref <- dimarConstructReferenceData(mtx)
#' sim <- dimarAssignPattern(ref, coef, mtx)

dimarAssignPattern <- function(ref, coef, mtx=NULL, npat=NULL) {
  if (is.null(npat)){
    if (dim(ref)[1]*dim(ref)[2]<50000) {
      npat <- 20
    } else if (dim(ref)[1]*dim(ref)[2]<100000){
      npat <- 10
    } else{
      npat <- 5
    }
  }

  X <- dimarConstructDesignMatrix(ref)
  pat <- array(NA,c(dim(ref),npat))

  for (i in 1:npat){
    # probability of MV by logistic regression
    yhat <- exp(X$X%*%coef)/(1+exp(X$X%*%coef))
    p <- matrix(yhat,nrow=dim(ref)[1])
    # binomial draw
    r <- matrix(runif(length(yhat)),nrow=dim(p)[1])
    ind <- which(r<p & !is.na(ref),arr.ind=T)
    # if ref has MV already, #MV = #MV of original mtx
    if (!is.null(mtx) & dim(ind)[1]>sum(is.na(mtx))) {
      ind <- ind[sample(1:length(ind),sum(is.na(mtx))),]
    }
    # Assign NA to ref
    pat1 <- ref
    pat1[ind] <- NA
    # if protein not measured at all, randomly assign one data point
    if (any(rowSums(is.na(pat1))==dim(pat1)[2])) {
      allna <- which(rowSums(is.na(pat1))==dim(pat1)[2])
      allna <- cbind(allna,sample(1:dim(pat1)[2],length(allna)))
      pat1[allna] <- ref[allna]
    }
    pat[,,i] <- pat1
  }
  print(paste(npat,'patterns of MVs are assigned.'))
  return(pat)
}
