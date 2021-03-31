dimar <- function(mtx) {
  
if (is.character(mtx)) {
  filename <- mtx
  mtx <- dimarReadInMaxQuant(mtx)
} else { filename <- {} }

coef <- dimarLearnPattern(mtx)
ref <- dimarConstructReferenceData(mtx)
sim <- dimarAssignPattern(ref, coef, mtx)

Imputations <- dimarDoImputations(sim, c('impSeqRob', 'ppca', 'imputePCA'))
Performance <- dimarEvaluatePerformance(Imputations, ref, sim, 'RMSE', TRUE)
Imp <- dimarDoOptimalImputation(mtx, rownames(Performance))

write.csv(Imp, file=paste0("Imp_", filename))
return(Imp)
}