#' dimarConstructRegularizationMatrix
#'
#' @description Generates regularization matrix: Expand the design matrix so each predictor has at least one positive entry, this exhibits the regression coefficients to get infinite
#' @return list of response vector, design matrix with regularization and predictor type
#' @param design list of response vector, design matrix and predictor type
#' @export dimarConstructRegularizationMatrix
#' @examples 
#' design <- dimarConstructDesignMatrix(mtx)
#' design <- dimarConstructRegularizationMatrix(design)

dimarConstructRegularizationMatrix <- function(design) {
  nrest <- sum(design$Xtype!=3 & design$Xtype!=2)
  Xreg <- matrix(0, 2*(length(design$Xtype)-nrest), length(design$Xtype))
  yreg <- matrix(0, 2*(length(design$Xtype)-nrest), 1)
  # median for the biological predictors
  Xreg[, 1:nrest] <- rep(median(design$X[,1:nrest]), dim(Xreg)[1])
  # 0/1 to regularize row/col coefficients
  id <- 1
  for (i in (nrest+1):dim(Xreg)[2]) {
    Xreg[id:(id+1), i] = 1
    yreg[id+1] = 1
    id = id+2
  }
  design$y <- c(design$y, yreg)
  design$X <- rbind(design$X, Xreg)
  return(design)
}
