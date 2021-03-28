#' dimarConstructReferenceData
#'
#' @description Constructs reference dataset
#' @return Reference dataset
#' @param mtx Quantitative matrix
#' @param cut Variable cut
#' @export dimarConstructReferenceData
#' @examples Sample example to demonstrate the function
dimarConstructReferenceData <- function(mtx, cut=0.2) {
  if (!is.numeric(cut) || cut<0 || cut>100){
    warning(paste('dimarConstructReferenceData.R: Variable cut =',cut,'is not supported. Check here. Used cut = 0.2 instead.'))
  }
  if (cut>1) {
    cut = cut/100
  }

  # sort by #NA
  mtx <- mtx[order(rowSums(is.na(mtx))),]
  nasum <- rowSums(is.na(mtx))
  # Assign mtx1 until cut
  nacut <- nasum[ceiling(dim(mtx)[1]*cut)]
  mtx1 <- mtx[nasum <= nacut,]
  idx1 <- 1:dim(mtx1)[1]
  idx2 <- (dim(mtx1)[1]+1):dim(mtx)[1]

  # Take mtx1 multiple times and normalize to protein means/std
  idxnew <- c()
  while (length(idxnew) < length(idx2)) {
    if (length(idxnew) + length(idx1) < length(idx2)) {
      idxnew = c(idxnew, idx1)
    } else {
      idxnew = c(idxnew, sample(length(idx1), length(idx2)-length(idxnew)))
    }
  }
  mtx1 <- rbind(mtx1, (mtx[idxnew,]-rowMeans(mtx[idxnew,],na.rm=T)) / rowSds(mtx[idxnew,], na.rm=T) * rowSds(mtx[idx2,], na.rm=T) + rowMeans(mtx[idx2,], na.rm=T) )
  mtx <- mtx1[order(rowSums(is.na(mtx1))),]
  print('Reference data is constructed.')
  return(mtx)
  }
