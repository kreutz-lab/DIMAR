#' dimarConstructDesignMatrix
#'
#' @description Constructs design matrix with the factorial row/col predictors and the mean intensity
#' @return list of response vector, design matrix and vector of predictor type
#' @param mtx Quantitative matrix
#' @param group a vector describing the assigment of each sample (column) of mtx to a group. By default mtx is split into two equally large groups
#' @param rowCoefByGroup by default TRUE. Thereby row coefficients are estimated for each group seperatly.
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

dimarConstructDesignMatrix <- function(mtx, group = rep(c(1,2), each = ncol(mtx)/2), rowCoefByGroup = T) {
  message("dimarConstructDesignMatrix by default generates row Coefficients for each group. By Default it is assumed the data contains two equally sizes groups, 
          with the first half of columns refering to samples of the first group. A costom vector can be passed to the function if another group design is needed.
          Missing Value Patterns may contain group information that are relevant if imputation is performed and are only kept if group specific row coefficients are used.")
  X <- matrix(0L, nrow = nrow(mtx)*ncol(mtx), ncol = nrow(mtx)*length(unique(group)) + ncol(mtx) + 2)
  # Intercept
  X[,1] <- rep(1, nrow(mtx)*ncol(mtx))
  Xname <- 'Intercept'
  Xtype <- 0
  # Mean intensity
  X[,2] <- scale(rep(rowMeans(mtx, na.rm = TRUE), ncol(mtx)))
  X[which(is.na(X[,2])),2]  <- min(X[,2], na.rm=T) ##needed if fitted by group. Some proteins have no values in one group but are to be kept. 
  Xname[2] <- 'mean'
  Xtype[2] <- 1
  #Assign row coefficients for each group. Only then group differences in the missing value pattern can be reflected in the coefficients
  if(rowCoefByGroup == F)
    group <- rep(1,ncol(mtx))
  row <- paste0(rep(rownames(mtx), ncol(mtx)), "#",rep(unique(group), each = ncol(mtx)/length(unique(group))*nrow(mtx)))
  col <- rep(1:ncol(mtx), each = nrow(mtx))
  
  for (i_col in 1:ncol(mtx)) {
    X[col == i_col, i_col + 2] <- 1
    Xname <- c(Xname, colnames(mtx)[i_col])
    Xtype <- c(Xtype, 2)
  }
  for (i_row in 1:(nrow(mtx)*length(unique(group)))) {
    rowID <- unique(row)[i_row]
    X[row == rowID, i_row + 2 + ncol(mtx)] <- 1  
    Xname <- c(Xname, rowID)  
    Xtype <- c(Xtype, 3)  
  }
  colnames(X) <- Xname
  
  y <- as.numeric(is.na(mtx))
  #return(list(y,X,Xtype))
  design <- list(y, X,Xtype)
  names(design) <- c('y', 'X','Xtype')
  return(design)
}
