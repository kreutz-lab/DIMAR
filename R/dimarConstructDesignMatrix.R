#' dimarConstructDesignMatrix
#'
#' @description Constructs design matrix
#' @return Design matrix
#' @param mtx Quantitative matrix
#' @export dimarConstructDesignMatrix
#' @examples Sample example to demonstrate the function
dimarConstructDesignMatrix <- function(mtx) {
  X <- matrix(0L, nrow=dim(mtx)[1]*dim(mtx)[2], ncol=dim(mtx)[1]+dim(mtx)[2]+2)
  # Intercept
  X[,1] <- rep(1, dim(mtx)[1]*dim(mtx)[2])
  Xname <- 'Intercept'
  Xtype <- 0
  # Mean intensity
  X[,2] <- scale(rep(rowMeans(mtx, na.rm=TRUE), dim(mtx)[2]))
  Xname[2] <- 'mean'
  Xtype[2] <- 1
  row <- rep(1:dim(mtx)[1], dim(mtx)[2])
  col <- rep(1:dim(mtx)[2], each=dim(mtx)[1])
  for (i in 1:dim(mtx)[2]) {
    X[col == i, i+2] <- 1
    Xname <- c(Xname, paste('Col',i))
    Xtype <- c(Xtype, 2)
  }
  for (i in 1:dim(mtx)[1]) {
    X[row==i, i+2+dim(mtx)[2]] <- 1
    Xname <- c(Xname, paste('Row', i))
    Xtype <- c(Xtype, 3)
  }
  y <- as.numeric(is.na(mtx))
  #return(list(y,X,Xtype))
  design <- list(y, X,Xtype)
  names(design) <- c('y', 'X','Xtype')
  return(design)
}
