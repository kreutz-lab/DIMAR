source("GetLib.R")
source("DoImputationsR.R")

DoOptimalImputation <- function(mtx,method,lib=NULL) {
  
m<-0
Imp <- NULL
while (is.null(Imp) || any(is.na(Imp))) {
  m<-m+1
  eval(parse(text=paste('require(',GetLib(method[m]),')')))
  Imp <- DoImputationsR(mtx,method[m],GetLib(method[m]))
}
Imp <- as.matrix(Imp[,1:dim(mtx)[2]])
print(paste('Imputation with',method[m],'done.'))
return(Imp)
}

