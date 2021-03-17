GetLib <- function(method) {
  
if (method=='midastouch' || method=='rf' || method=='mean' || method=='norm' || method=='ri' || method=='pmm' || method=='sample' || method=='cart')
  { lib <-'mice' }
if (method=='knn' || method=='impute.knn') 
  { lib <-'impute'}
if (method=='impnorm' || method=='imp.norm')
  { lib <-'norm'}
if (method=='Amelia' || method=='amelia')
  {lib = 'Amelia'}
if (method=='regression' || method=='aregImpute')
  { lib <-'Hmisc'}
if (method=='ppca' || method=='bpca' || method=='nipals' || method=='nlpca' || method=='svd' || method=='svdImpute')
  { lib <-'pcaMethods'}
if (method=='MinDet' || method=='KNN' || method=='MinProb' || method=='QRILC')
  { lib <-'imputeLCMD'}
if (method=='SVTApproxImpute' || method=='SVTImpute' || method=='SVDImpute' || method=='kNNImpute' || method=='lmImpute')
  { lib <-'imputation'}
if (method=='missForest')
  { lib <-'missForest'}
if (method=='softImpute') { lib <-'softImpute'}
if (method=='irmi') { lib <-'VIM'}
if (method=='Seq' || method=='SeqRob' || method=='impSeq' || method=='impSeqRob') 
  { lib <-'rrcovNA'}
if (method=='MIPCA' || method=='imputePCA')
  { lib <-'missMDA'}
if (method=='mi') { lib <-'mi'}
if (method=='knnImputation'){ lib <-'DMwR'}
if (method=='GMSimpute' || method=='GMSLasso')
  { lib <-'GMSimpute'}
if (!exists('lib')) {
  warning(paste('GetLib.R:',method,'is not implemented in DIMA and ignored for imputation. Check your spelling, or add respective Rcode to DoImputationsR.R and its library to GetLib.R.'))
}
return(lib)
}