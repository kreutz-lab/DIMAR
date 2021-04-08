#' dimar
#' @description DIMAR: Data-driven selection of an imputation algorithm in R
#' @param mtx Data matrix or input file (.txt, .xls, .xlsx, .csv or .tsv)
#' @param pattern Search pattern for specifying sample names read in as default data, if not specified the user will be asked
#' @export dimar dimar
#' @examples 
#' Imp <- dimar(mtx)
#' Imp <- dimar(file)
#' Imp <- dimar('proteinGroups.txt','LFQ')

dimar <- function(mtx,pattern="^Int") {
  
if (is.character(mtx)) {
  file <- mtx
  ext <- strsplit(basename(file), split="\\.")[[1]][-1]
  mtx <- read.table(file,header=T,sep="\t", allowEscapes = TRUE, check.names = FALSE)
} else { filename <- {} }
if (class(mtx) == 'SingleCellExperiment' || class(mtx) == 'SummarizedExperiment') {
  mtx <- as.matrix(assay(mtx))
}
if (is.character(pattern)) {
  mtx <- as.matrix(mtx[, grepl(pattern, names(mtx))])
  print(paste("Data is reduced to columns which include",pattern,"in their sample names."))
}

mtx <- dimarMatrixPreparation(mtx)

coef <- dimarLearnPattern(mtx)
ref <- dimarConstructReferenceData(mtx)
sim <- dimarAssignPattern(ref, coef, mtx)

Imputations <- dimarDoImputations(sim, c('impSeqRob', 'ppca', 'imputePCA'))
Performance <- dimarEvaluatePerformance(Imputations, ref, sim, 'RMSE', TRUE)
Imp <- dimarDoOptimalImputation(mtx, rownames(Performance))

write.table(Imp, file=file.path(dirname(file),paste0("Imp_",basename(file))),sep="\t")
print(paste('Imputation is written in',file.path(dirname(file),paste0("Imp_",basename(file)))))
return(Imp)
}