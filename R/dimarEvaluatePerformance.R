#' dimarEvaluatePerformance
#'
#' @description Evaluates performance of imputation algorithms.
#' @return Data frame containing the following performance measures for
#' each imputation method: Deviation, RMSE, RSR, p-Value_F-test, Accuracy, PCC, and in case
#' of RMSEttest=TRUE the RMSE t-test result
#' @param Imputations Imputed data set(s)
#' @param ref Reference data
#' @param sim Simulated patterns of MVs
#' @param rankby Performance measure which should serve as rank criterion
#' @param RMSEttest flag if RMSE of ttest should be calculated
#' @param group indices for ttest
#' @export dimarEvaluatePerformance dimarEvaluatePerformance
#' @examples
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' coef <- dimarLearnPattern(mtx)
#' ref <- dimarConstructReferenceData(mtx)
#' sim <- dimarAssignPattern(ref, coef, mtx)
#' Imputations <- dimarDoImputations(sim, c('impSeqRob', 'ppca', 'imputePCA'))
#' Performance <- dimarEvaluatePerformance(Imputations, ref, sim)

dimarEvaluatePerformance <- function(Imputations, ref, sim, rankby = 'RMSE',
                                     RMSEttest = TRUE, group = 'cluster') {
  # Initialize performance arrays
  Dev <- array(NA, c(dim(Imputations[[1]])[3], length(Imputations)))
  RMSE <- Dev
  RSR <- Dev
  pF <- Dev
  Acc <- Dev
  PCC <- Dev
  rank <- Dev
  if (RMSEttest) {
  RMSEt <- Dev
  ttest <- NULL
  ttesti <- NULL
  if (group[1] == 'cluster' | is.null(group)) {
    h <- amap::hcluster(t(ref))
    group <- stats::cutree(h, k = 2)
    if (sum(group == 1) < 2 | sum(group == 2) < 2) { RMSEttest <- FALSE}
  }
  }

  for (p in 1:dim(Imputations[[1]])[3]) { # loop over #patterns
    if (RMSEttest) {
      for (t in 1:dim(ref)[1]) {
        htest = stats::t.test(ref[t, group == 1], ref[t, group == 2])
        ttest[t] <- htest$statistic
      }
      ttest[!is.finite(ttest)] <- NULL
    }
    for (a in 1:length(Imputations)) { # loop over imputation algorithms
      im <- Imputations[[a]][, , p]
      if (!any(is.na(im))) {
          ndata = sum(is.na(sim[, , p]) & !is.na(ref))
          Diff <- im - ref
          Dev[p,a] <- sum(abs(Diff), na.rm = TRUE)/ndata
          RMSE[p,a] <- sqrt(sum(Diff^2, na.rm = TRUE)/ndata)
          RSR[p,a] <- RMSE[p,a]/stats::sd(ref, na.rm = TRUE)
          pF[p,a] <- stats::var.test(im,ref)$p.value
          Acc[p,a] <- length(which(abs(Diff/ref) < 0.05))/dim(ref)[1]/dim(ref)[2]*100
          PCC[p,a] <- stats::cor(as.vector(im), as.vector(ref))
          if (RMSEttest) {
            for (t in 1:dim(ref)[1]) {
              htesti <- stats::t.test(im[t,group == 1], im[t,group == 2])
              ttesti[t] <- htesti$statistic
            }
            ttesti[!is.finite(ttesti)] <- NULL
            RMSEt[p,a] <- sqrt(sum((ttest - ttesti)^2, na.rm = TRUE) / max(length(ttest), length(ttesti)))
          }
      }
    }
    rank[p,] <- eval(parse(text = paste('order(', rankby, '[p,])', sep = "")))
  }
  if (p == 1) {
  if (RMSEttest) {
    Performance <- data.frame(Dev, RMSE, RSR, pF, Acc, PCC, RMSEt)
  } else {
    Performance <- data.frame(Dev, RMSE, RSR, pF, Acc, PCC)
  }
  } else {
  rank <- colMeans(rank)
  if (RMSEttest) {
    Performance <- data.frame(colMeans(Dev), colMeans(RMSE), colMeans(RSR), colMeans(pF),
                              colMeans(Acc), colMeans(PCC), colMeans(RMSEt))
  } else {
    Performance <- data.frame(colMeans(Dev), colMeans(RMSE), colMeans(RSR), colMeans(pF),
                              colMeans(Acc), colMeans(PCC))
  }
  }
  Performance <- Performance[rank,]
  rownames(Performance) <- names(Imputations)
  if (RMSEttest) {
  colnames(Performance) <- c('Deviation', 'RMSE', 'RSR', 'p-Value_F-test', 'Accuracy',
                             'PCC', 'RMSEttest')
  } else {
  colnames(Performance) <- c('Deviation', 'RMSE', 'RSR', 'p-Value_F-test', 'Accuracy',
                             'PCC')
  }
  return(Performance)
}
