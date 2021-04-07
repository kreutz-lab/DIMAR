#' dimar
#' @description DIMAR: Data-driven selection of an imputation algorithm in R
#' @param mtx Data matrix or input file (.txt, .xls, .xlsx, .csv or .tsv)
#' @param pattern Search pattern for specifying sample names read in as default data, if not specified the user will be asked
#' @export dimar dimar
#' @examples 
#' Imp <- dimar(mtx)
#' Imp <- dimar(file)
#' Imp <- dimar('proteinGroups.txt','LFQ')

dimar <- function(mtx,pattern=NULL) {
  
if (is.character(mtx)) {
  file <- mtx
  ext <- strsplit(basename(file), split="\\.")[[1]][-1]
  if (ext=='tsv'){
    mtx <- read.table(file='metabolomics.tsv',header=T)
  } else if (ext=='csv'){
    mtx <- read.csv(filename, allowEscapes = TRUE, check.names = FALSE,sep = "\t")
  } else if (ext=='txt'){
    mtx <- read.delim(file)
  } else {
    warning('Input file extension unknown. Expand code here or check input.')
  }
} else { filename <- {} }
if (class(mtx) == 'SingleCellExperiment' || class(mtx) == 'SummarizedExperiment') {
  mtx <- as.matrix(assay(mtx))
}
if (is.character(pattern)) {
  mtx <- as.matrix(mtx[, grepl(pattern, names(mtx))])
}

mtx <- dimarMatrixPreparation(mtx)

coef <- dimarLearnPattern(mtx)
ref <- dimarConstructReferenceData(mtx)
sim <- dimarAssignPattern(ref, coef, mtx)

Imputations <- dimarDoImputations(sim, c('impSeqRob', 'ppca', 'imputePCA'))
Performance <- dimarEvaluatePerformance(Imputations, ref, sim, 'RMSE', TRUE)
Imp <- dimarDoOptimalImputation(mtx, rownames(Performance))

write.csv(Imp, file=paste0("Imp_", file))
return(Imp)
}