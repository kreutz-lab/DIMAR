
EvaluatePerformance <- function(Imputations,ref,sim,rankby='RMSE',RMSEttest='false',group='cluster') {
  
# Initialize performance arrays
Dev <- array(NA,c(dim(Imputations$Imp)[3],dim(Imputations$Imp)[4]))
RMSE <- Dev
RSR <- Dev
pF <- Dev
Acc <- Dev
PCC <- Dev

for (p in 1:dim(Imputations$Imp)[3]) { # loop over #patterns
  for (a in 1:dim(Imputations$Imp)[4]) { # loop over imputation algorithms
    im <- Imputations$Imp[,,p,a]
    if (!any(is.na(im))) {
        ndata = sum(is.na(sim[,,p]) & !is.na(ref))
        Diff <- im-ref
        Dev[p,a] <- sum(abs(Diff),na.rm=T)/ndata
        RMSE[p,a] <- sqrt(sum(Diff^2,na.rm=T)/ndata)
        RSR[p,a] <- RMSE[p,a]/sd(ref,na.rm=T)
        pF[p,a] <- var.test(im,ref)$p.value
        Acc[p,a] <- length(which(abs(Diff/ref)<0.05))/dim(ref)[1]/dim(ref)[2]*100
        PCC[p,a] <- cor(as.vector(im),as.vector(ref))
    }
  }
}

}