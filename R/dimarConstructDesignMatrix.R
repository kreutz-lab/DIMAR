#' dimarConstructDesignMatrix
#'
#' @description Constructs design matrix with the factorial row/col predictors and the mean intensity
#' @return list of response vector, design matrix and vector of predictor type
#' @param mtx Quantitative matrix
#' @param ind indices of rowSubset
#' @param group a vector describing the assigment of each sample (column) of mtx to a group. 
#' By default mtx is calculated using hierarchical clustering.
#' Missing Value Patterns may contain group information that are relevant if imputation is performed and are only kept if group specific row coefficients are used.
#' @param rowCoefByGroup by default TRUE. Thereby row coefficients are estimated for each group seperatly.
#' @param orderCoefByName Keep names of original input matrix in estimated coefficients
#' @param DE_idx If groundtruth is known, give this parameter to ensure DE proteins from input also receive row coefficients of DE proteins
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

dimarConstructDesignMatrix <- function(mtx, ind = 1:nrow(mtx), 
                                       group =  'cluster', 
                                       rowCoefByGroup = T, 
                                       orderCoefByName = F, 
                                       DE_idx = NULL) {
  if(rowCoefByGroup ==T){
  if(is.null(group))
    group <- 'cluster'
  if (group[1] == 'cluster') {
    message("By default DIMAR is using hierarchical clustering (amap::hcluster())
            to assign samples to two groups. If group assigments are known,
            provide a vector of group assigments to be passed to the 
            dimarConstructDesignMatrix function.")
    h <- amap::hcluster(t(mtx))
    group <- stats::cutree(h, k = 2)
  }
  }else{ 
    message("Row coefficients are estimated for all groups together. 
            This may result in averaging out the amount of missing values over all groups.
            If you want to estimate row coefficients for each group seperatly, set rowCoefByGroup to TRUE.")
    group <- rep(1,ncol(mtx))}
  
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
  

  rowIDs <- paste0("row",ind)
  colIDs <- paste0("col",1:ncol(mtx))
  #make rowIDs. If groundtruth is known, information about grundtruth of a protein can be stored in it's row coefficient name
  if(!is.null(DE_idx)){
    DE_idx_chunk <- which(ind %in% DE_idx)
    rowIDs[DE_idx_chunk] <- paste0(rowIDs[DE_idx_chunk],"_DE")
  }
  #if you want to use the original names from mtx to keep the exact order
  if(orderCoefByName){
    if(is.null(rownames(mtx)))
      stop("The input matrix has no rownames.
      Please provide rownames to keep the order of the coefficients,
           or set orderCoefByName to FALSE.")
    rowIDs <- rownames(mtx)
    if(is.null(colnames(mtx)))
       stop("The input matrix has no colnames.
      Please provide rownames to keep the order of the coefficients,
           or set orderCoefByName to FALSE.")
    colIDs <- colnames(mtx)
  }
  row <- paste0(rep(rowIDs, ncol(mtx)), "#",rep(unique(group), each = ncol(mtx)/length(unique(group))*nrow(mtx)))
  col <- rep(colIDs, each = nrow(mtx))

  for (i_col in 1:ncol(mtx)) {
    colID <- unique(col)[i_col]
    X[col == colID, i_col + 2] <- 1
    Xname <- c(Xname, colID)
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
