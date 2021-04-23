#' dimarDoImputationsR
#'
#' @description Helper function for DoImputations. Actually applies the imputation functions
#' of various R packages to a quantitative matrix.
#' @return Quantitaive matrix, on which imputation is performed
#' @param mtx Quantitative matrix
#' @param method Imputation method
#' @param lib R package of imputation method
#' @export dimarDoImputationsR
#' @examples
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' I <- dimarDoImputationsR(mtx,'impSeq','rrcovNA')

dimarDoImputationsR <- function(mtx, method = NULL, lib = NULL) {
  if (is.null(method)) {
    warning('DoImputations.R: No method given. No imputation performed.')
    Imp <- {}
  }
  if (is.null(lib)) {
    lib <- dimarGetLib(method)
  }
  tryCatch({
    if (lib == 'pcaMethods') {
      mtx[is.na(mtx)] <- NA
      if (sum(rowSums(is.na(mtx)) >= ncol(mtx)) > 0) {
        Imp <- {}
      } else {
        I <- pcaMethods::pca(mtx, method = method)
        Imp <- pcaMethods::completeObs(I)
      }
    } else if (lib == 'impute') {
      I <- impute::impute.knn(as.matrix(mtx))
      Imp <- I$data
    } else if (lib == 'norm') {
      s <- norm::prelim.norm(mtx)
      thetahat <- norm::em.norm(s)
      norm::rngseed(1)
      Imp <- norm::imp.norm(s, thetahat, mtx)
    } else if (lib == 'missMDA') {
      f <- get(method)
      imp <- f(data.frame(mtx), nboot = 1)
      if (method == 'MIPCA') {
        Imp <- imp$res.imputePCA
      } else {
        Imp <- imp$completeObs
      }
    } else if (lib == 'rrcovNA') {
      if (grepl('imp',method)) {
        f <- get(method)
      } else {
        f <- get(paste('imp', method, sep = ''))
      }
      Imp <- f(mtx)
      if (grepl('SeqRob',method)) {
        Imp <- Imp$x
      }
    } else if (lib == 'VIM') {
      f <- get(method)
      Imp <- f(as.matrix(mtx))
    } else if (lib == 'softImpute') {
      f <- softImpute::softImpute(as.matrix(mtx))
      Imp <- softImpute::complete(as.matrix(mtx),f)
    } else if (lib == 'imputeLCMD') {
      if (method == 'QRILC') {
        Imp <- imputeLCMD::impute.QRILC(as.matrix(mtx))[[1]]
      } else {
        f <- get(paste('impute.', method, sep = ''))
        Imp <- f(as.matrix(mtx))
      }
    } else if (lib == 'imputation') {
      f <- get(method)
      if (grepl('SVD', method)) {
        Imp <- f(mtx, k = 3)$x
      } else if (grepl('kNN', method)) {
        Imp <- f(mtx, k = 3)$x
      } else if (grepl('SVT', method)) {
        Imp <- f(mtx, lambda = 3)$x
      } else {
        I <- f(mtx)
        Imp <- I$x
      }
    } else if (lib == 'mice') {
      colnames(mtx) <- c(paste0("X",1:dim(mtx)[2]))
      I <- mice::mice(mtx, m = 1, method = method)
      Imp <- mice::complete(I)
    } else if (lib == 'Amelia') {
      # if isSymmetric, R aborts without error message
      if (isSymmetric(mtx)) {
        f <- R.utils::withTimeout({Amelia::amelia(mtx, m = 1)}, timeout = 1, cpu = 100, elapsed = 3600) # of all pride mtx the max time of amelia was cpu=1
        Imp <- f$imputations[[1]]
      } else {
        Imp <- NULL
      }
    } else if (lib == 'missForest') {
      f <- missForest::missForest(mtx)
      Imp <- f$ximp
    } else if (lib == 'Hmisc') {
      colnames(mtx) <- c(paste0("X", 1:dim(mtx)[2]))
      formula <- '~ X1'
      for (j in 2:dim(mtx)[2]) {
        formula <- paste(formula,' +X', j, sep = '')
      }
      if (method == 'aregImpute') {
        f <- Hmisc::aregImpute(stats::as.formula(formula), data = data.frame(mtx), n.impute = 1, type = "pmm")
      } else {
        f <- Hmisc::aregImpute(stats::as.formula(formula), data = data.frame(mtx), n.impute = 1, type = method)
      }
      Imp <- array(unlist( Hmisc::impute.transcan(f, imputation = TRUE, data = data.frame(mtx),
                                                  list.out = TRUE)), dim = dim(mtx ))

    } else if (lib == 'DMwR') {
      Imp <- DMwR::knnImputation(mtx)
    } else if (lib == 'mi') {
      I <- mi::mi(data.frame(mtx), n.chains = 1)
      Imp <- mi::complete(I)[1:length(mtx)]
    } else if (lib == 'GMSimpute') {
      Imp <- GMSimpute::GMS.Lasso(mtx, log.scale = TRUE, TS.Lasso = TRUE)
    } else {
      Imp <- NULL
      stop(paste('dimarDoImputationsR.m: library', lib, 'is not recognized. Expand code here.'))
    }
  }, error = function(e) {
    Imp <<- NULL
    warning(paste('dimarDoImputationsR.R: Error in R package', lib, 'within algorithm', method,':', conditionMessage(e)))
  })
  return(Imp)
}


