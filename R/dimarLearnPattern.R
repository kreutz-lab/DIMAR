#' dimarLearnPattern
#'
#' @description Learns missing value occurence pattern from quantitative matrix by a logistic regression model
#' @return Logistic regression coefficients describing missing value pattern
#' @param mtx Quantitative matrix
#' @export dimarLearnPattern
#' @examples
#' #' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' coef <- dimarLearnPattern(mtx)

dimarLearnPattern <- function(mtx) {
  # Subsample indices
  if (nrow(mtx) > 1000) {
    nsub <- ceiling(nrow(mtx) / 1000)
    indrand <- sample(1:nrow(mtx), nrow(mtx))
    npersub <- ceiling(nrow(mtx) / nsub)
  } else {
    nsub <- 1
    ind <- 1:nrow(mtx)
  }

  for (i in 1:nsub) {
    #  i <- 1
    if (nsub > 1) {
      ind <- indrand[(npersub*(i - 1) + 1):(npersub*i)]
    }
    design <- dimarConstructDesignMatrix(mtx[ind,])
    design <- dimarConstructRegularizationMatrix(design)

    #fit <- stats::glm.fit(X,y,family=stats::binomial(),weights=rep(1,dim(X)[1]))
    fit <- stats::glm.fit(design$X, design$y, family = stats::binomial())
    if (i == 1) {
      coef <- stats::coefficients(fit)
    } else {
      coef <- rbind(coef, stats::coefficients(fit))
    }
  }
  # sort row coefficients, intensity/column coefficients are set to mean over nsub (for loop)
  if (nsub > 1) {
    coef <- c(colMeans(coef[, design$Xtype != 3]), sort(coef[, design$Xtype == 3]))
  }

  print('Pattern of MVs is learned by logistic regression.')
  return(coef)
}
