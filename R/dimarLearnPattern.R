#' dimarLearnPattern
#'
#' @description Learns missing value occurence pattern from quantitative matrix by a logistic regression model
#' @return Logistic regression coefficients describing missing value pattern
#' @param mtx Quantitative matrix
#' @param orderCoefByName default: FALSE. If coefficients should be ordered exactly as in input mtx set TRUE. Then, eg row 100 of a new matrix will 
#' be simulated with the row coefficient of row 100 of the input matrix
#' @param DE_idx if groundtruth is known, give the row index of differentially expressed proteins
#' @export dimarLearnPattern
#' @examples
#' #' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' coef <- dimarLearnPattern(mtx)

dimarLearnPattern <- function(mtx, orderCoefByName = F, DE_idx = NULL) {
  # Subsample indices
  if (nrow(mtx) > 500) {
    nsub <- ceiling(nrow(mtx) / 500)
    indrand <- sample(1:nrow(mtx), nrow(mtx))
    npersub <- ceiling(nrow(mtx) / nsub)
    coef.lst <- list() #intialise as list in order to keep names of rows
  } else {
    nsub <- 1
    ind <- 1:nrow(mtx)
  }

  for (i in 1:nsub) {
    #  i <- 1
    if (nsub > 1) {
      if (i==nsub) {
        ind <- indrand[(npersub*(i-1)+1):length(indrand)]
      } else {
        ind <- indrand[(npersub*(i - 1) + 1):(npersub*i)]
      }
    }
    design <- dimarConstructDesignMatrix(mtx[ind,], ind, DE_idx = DE_idx, orderCoefByName = orderCoefByName)
    design <- dimarConstructRegularizationMatrix(design)

    #fit <- stats::glm.fit(X,y,family=stats::binomial(),weights=rep(1,dim(X)[1]))
    fit <- stats::glm.fit(design$X, design$y, family = stats::binomial())
    if (nsub==1){
      coef <- stats::coefficients(fit)
    } else {
      #coef[i,1:length(coefficients(fit))] <- stats::coefficients(fit)
      coef.lst[[i]] <- stats::coefficients(fit)
    }
  }
  # sort row coefficients, intensity/column coefficients are set to mean over nsub (for loop)
  if (nsub > 1) {
    idx <- which(design$Xtype!=3)
    coef.mtx <- dplyr::bind_rows(coef.lst)
    notRowCoef <- colMeans(coef.mtx[,idx])
    notRowCoef["Intercept"] <- 0 #intercept is substracted from row id coefficients
    rowCoef.Intercept <- apply(coef.mtx[,-idx],2,
                               function(x) as.numeric(x + unlist(coef.mtx[,1])))
    rowCoef <- colMeans(rowCoef.Intercept, na.rm=T)
    
    rowCoef <- reorderRowCoef(rowCoef,mtx, orderCoefByName = orderCoefByName)
    
    coef <- c(notRowCoef,rowCoef)
  }

  print('Pattern of MVs is learned by logistic regression.')
  return(coef)
}

reorderRowCoef <- function(rowCoef,mtx, orderCoefByName = F){
  rowGroupVector <- gsub(".*#", "", names(rowCoef))
  rowCoef_ordered <- c()
  for (g in unique(rowGroupVector)){
    rowCoef.group <- rowCoef[grep(paste0("#",g), names(rowCoef))]
    protNames.group <- gsub("#.*", "", names(rowCoef.group))
    if(orderCoefByName)
      rowCoef.group <- rowCoef.group[order(match(protNames.group,rownames(mtx)))] #reorder row Coefficients so they resemble order in mtx
    rowCoef_ordered <- c(rowCoef_ordered, rowCoef.group)
  }
  return(rowCoef_ordered)
}