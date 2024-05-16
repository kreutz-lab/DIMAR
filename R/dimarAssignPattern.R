#' dimarAssignPattern
#'
#' @description Assigns patterns of MVs to reference data set with a logistic regression model
#' @return Reference dataset with applied pattern
#' @param ref Reference dataset
#' @param coef Logistic regression coefficients of missing value pattern
#' @param mtx Quantitative matrix
#' @param npat Number of patterns
#' @param group a vector describing the assigment of each sample (column) of mtx to a group.
#' if not known, group is calculated using hierarchical clustering
#' @param dontAllowOnlyNA Default: TRUE ; Proteins with only NA will be assigned a random value.
#' If DIMAR is used to simulate by group, set this to FALSE. Proteins with only missing values in one of the groups are common and important for a realistic missing value pattern.
#' @export dimarAssignPattern
#' @examples
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' coef <- dimarLearnPattern(mtx)
#' ref <- dimarConstructReferenceData(mtx)
#' sim <- dimarAssignPattern(ref, coef, mtx)

dimarAssignPattern <- function(ref, coef, mtx = NULL, npat = NULL, group = NULL, dontAllowOnlyNA = TRUE) {

  if (is.null(npat)) {
    if (nrow(ref)*ncol(ref) < 50000) {
      npat <- 20
    } else if (nrow(ref)*ncol(ref) < 100000) {
      npat <- 10
    } else{
      npat <- 5
    }
  }

  X <- dimarConstructDesignMatrix(ref, group = group)
  pat <- array(NA,c(dim(ref),npat))

  for (i in 1:npat) {
    # probability of MV by logistic regression
    yhat <- exp(X$X %*% coef)/(1 + exp(X$X %*% coef))
    p <- matrix(yhat,nrow = nrow(ref))
    # binomial draw
    r <- matrix(stats::runif(length(yhat)),nrow = dim(p)[1])
    ind <- which(r < p & !is.na(ref),arr.ind = TRUE)
    # if ref has MV already, #MV = #MV of original mtx
    if (!is.null(mtx) & nrow(ind) > sum(is.na(mtx))) {
      ind <- ind[sample(1:nrow(ind),sum(is.na(mtx))),]
    }
    # Assign NA to ref
    pat1 <- ref
    pat1[ind] <- NA
    if(dontAllowOnlyNA){
    # if protein not measured at all, randomly assign one data point
      if (any(rowSums(is.na(pat1)) == ncol(pat1))) {
        allna <- which(rowSums(is.na(pat1)) == ncol(pat1))
        allna <- cbind(allna,sample(1:ncol(pat1),length(allna),replace=TRUE))
        pat1[allna] <- ref[allna]
      }
    }
    pat[,,i] <- pat1
  }
  message(paste(npat,'patterns of MVs are assigned.'))
  return(pat)


  }

