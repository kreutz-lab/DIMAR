#' dimarConstructDesignMatrix
#'
#' @description Constructs design matrix with the factorial row/col predictors and the mean intensity
#' @return list of response vector, design matrix and vector of predictor type
#' @param mtx Quantitative matrix
#' @export dimarConstructDesignMatrix
#' @examples
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' design <- dimarConstructDesignMatrix(mtx)
#' design <- dimarConstructRegularizationMatrix(design)
#' fit <- glm.fit(design$X, design$y, family=binomial())
#'
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' design <- dimarConstructDesignMatrix(mtx)
#' fit <- glm.fit(design$X, design$y, family=binomial())
#' yhat <- exp(design$X %*% coefficients(fit))/(1+exp(design$X %*% coefficients(fit)))

dimarConstructDesignMatrix <- function(mtx) {
  X <- matrix(0L, nrow = nrow(mtx)*ncol(mtx), ncol = nrow(mtx) + ncol(mtx) + 2)
  # Intercept
  X[,1] <- rep(1, nrow(mtx)*ncol(mtx))
  Xname <- 'Intercept'
  Xtype <- 0
  # Mean intensity
  X[,2] <- scale(rep(rowMeans(mtx, na.rm = TRUE), ncol(mtx)))
  Xname[2] <- 'mean'
  Xtype[2] <- 1
  row <- rep(1:nrow(mtx), ncol(mtx))
  col <- rep(1:ncol(mtx), each = nrow(mtx))
  for (i in 1:ncol(mtx)) {
    X[col == i, i + 2] <- 1
    Xname <- c(Xname, paste('Col',i))
    Xtype <- c(Xtype, 2)
  }
  for (i in 1:nrow(mtx)) {
    X[row == i, i + 2 + ncol(mtx)] <- 1
    Xname <- c(Xname, paste('Row', i))
    Xtype <- c(Xtype, 3)
  }
  y <- as.numeric(is.na(mtx))
  #return(list(y,X,Xtype))
  design <- list(y, X,Xtype)
  names(design) <- c('y', 'X','Xtype')
  return(design)
}
