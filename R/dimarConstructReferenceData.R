#' dimarConstructReferenceData
#'
#' @description Constructs reference dataset: Take data with least #MVs, enlarge data by proteins with least #MVs and normalize to all proteins
#' @return Reference dataset
#' @param mtx Quantitative matrix
#' @param cut minimum number of proteins kept for the reference data
#' @export dimarConstructReferenceData
#' @examples
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' coef <- dimarLearnPattern(mtx)
#' ref <- dimarConstructReferenceData(mtx)
#' sim <- dimarAssignPattern(ref, coef, mtx)

dimarConstructReferenceData <- function(mtx, cut = 0.2) {
  if (!is.numeric(cut) || cut < 0 || cut > 100) {
    warning(paste('dimarConstructReferenceData.R: Variable cut =', cut, 'is not supported. Check here. Used cut = 0.2 instead.'))
  }
  if (cut > 1) {
    cut = cut/100
  }

  # sort by #NA
  mtx <- mtx[order(rowSums(is.na(mtx))),]
  nasum <- rowSums(is.na(mtx))
  # Assign mtx1 until cut
  nacut <- nasum[ceiling(nrow(mtx)*cut)]
  mtx1 <- mtx[nasum <= nacut,]
  idx1 <- 1:nrow(mtx1)
  idx2 <- (nrow(mtx1) + 1):nrow(mtx)

  # Take mtx1 multiple times and normalize to protein means/std
  idxnew <- c()
  while (length(idxnew) < length(idx2)) {
    if (length(idxnew) + length(idx1) < length(idx2)) {
      idxnew = c(idxnew, idx1)
    } else {
      idxnew = c(idxnew, sample(length(idx1), length(idx2) - length(idxnew)))
    }
  }
  mtx1 <- rbind(mtx1, (mtx[idxnew,] - rowMeans(mtx[idxnew,], na.rm = TRUE)) /
                  matrixStats::rowSds(mtx[idxnew,], na.rm = TRUE) *
                  matrixStats::rowSds(mtx[idx2,], na.rm = TRUE) +
                  rowMeans(mtx[idx2,], na.rm = TRUE) )
  mtx <- mtx1[order(rowSums(is.na(mtx1))),]
  print('Reference data is constructed.')
  return(mtx)
  }
