dimar <- function(mtx,pattern<-NULL) {
  
if (is.character(mtx)) {
  file <- mtx
  ext <- strsplit(basename(file), split="\\.")[[1]][-1]
  if (ext=='tsv'){
    mtx <- read.table(file='metabolomics.tsv',header=T)
  } else if (ext=='csv'){
    mtx <- read.csv(filename, allowEscapes = TRUE, check.names = FALSE,sep = "\t")
  }
} else { filename <- {} }
if (class(mtx) == 'SingleCellExperiment' || class(mtx) == 'SummarizedExperiment') {
  mtx <- as.matrix(assay(mtx))
}
if (is.character(pattern)) {
  mtx <- as.matrix(mtx[, grepl(pattern, names(mtx))])
}

mtx <- dimarFiltering(mtx)

coef <- dimarLearnPattern(mtx)
ref <- dimarConstructReferenceData(mtx)
sim <- dimarAssignPattern(ref, coef, mtx)

Imputations <- dimarDoImputations(sim, c('impSeqRob', 'ppca', 'imputePCA'))
Performance <- dimarEvaluatePerformance(Imputations, ref, sim, 'RMSE', TRUE)
Imp <- dimarDoOptimalImputation(mtx, rownames(Performance))

write.csv(Imp, file=paste0("Imp_", file))
return(Imp)
}