#' dimar
#' @description DIMAR: Data-driven selection of an imputation algorithm in R
#' @param mtx Data matrix or input file (.txt, .xls, .xlsx, .csv or .tsv)
#' @param pattern Search pattern for specifying sample names read in as default data, if not specified the user will be asked
#' @param group vector of group indices for ttest (group==1 vs group==2) ['cluster'] as a default, clustering with 2 cluster is performed
#' @export dimar dimar
#' @examples 
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' Imp <- dimar(mtx)
#'
#' Imp <- dimar('proteinGroups.txt','LFQ')
#' 
#' Imp <- dimar('proteinGroups_PXD008893.txt','Intensity',c('PKB','PKC')))

dimar <- function(mtx,pattern="^Int",group='cluster') {

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
if (!group[1]=='cluster'){
  groupidx <- rep(0L,ncol(mtx))
  groupidx[grepl(group[1],colnames(mtx))]<-1
  groupidx[grepl(group[2],colnames(mtx))]<-2
  group <- groupidx
}
mtx <- dimarMatrixPreparation(mtx)

coef <- dimarLearnPattern(mtx)
ref <- dimarConstructReferenceData(mtx)
sim <- dimarAssignPattern(ref, coef, mtx)

Imputations <- dimarDoImputations(sim, c('impSeqRob', 'ppca', 'imputePCA'))
Performance <- dimarEvaluatePerformance(Imputations, ref, sim, 'RMSE', TRUE, group)
Imp <- dimarDoOptimalImputation(mtx, rownames(Performance))

if (exist("file")) {  # if file name is given, the imputed matrix is written into file
  write.table(Imp, file=file.path(dirname(file),paste0("Imp_",basename(file))),sep="\t")
  print(paste('Imputation is written in',file.path(dirname(file),paste0("Imp_",basename(file)))))
}
return(Imp)
}