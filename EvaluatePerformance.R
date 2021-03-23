
EvaluatePerformance <- function(Imputations,ref,sim,rankby='RMSE',RMSEttest=TRUE,group='cluster') {
  
# Initialize performance arrays
Dev <- array(NA,c(dim(Imputations$Imp)[3],dim(Imputations$Imp)[4]))
RMSE <- Dev
RSR <- Dev
pF <- Dev
Acc <- Dev
PCC <- Dev
rank <- Dev
if (RMSEttest) { 
  RMSEt <- Dev
  ttest<-NULL
  ttesti <- NULL 
  if (group=='cluster' | is.null(group)) {
    require(amap)
    h <- hcluster(t(ref))
    group <- cutree(h,k=2)
  }
}

for (p in 1:dim(Imputations$Imp)[3]) { # loop over #patterns
  if (RMSEttest) {
    for (t in 1:dim(Imputations$Imp)[1]) {
      htest = t.test(ref[t,group==1],ref[t,group==2])
      ttest[t] <- htest$statistic
    }
    ttest[!is.finite(ttest)] <- NULL
  }
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
        if (RMSEttest){
          for (t in 1:dim(ref)[1]) {
            htesti <- t.test(im[t,group==1],im[t,group==2])
            ttesti[t] <- htesti$statistic
          }
          ttesti[!is.finite(ttesti)] <- NULL
          RMSEt[p,a] <- sqrt(sum((ttest-ttesti)^2,na.rm=T) /max(length(ttest),length(ttesti)))
        }
    }
  }
  rank[p,] <- eval(parse(text=paste('order(',rankby,'[p,])',sep="")))
}
if (p==1){
  if (RMSEttest){
    Performance <- data.frame(Dev,RMSE,RSR,pF,Acc,PCC,RMSEt)
  } else {
    Performance <- data.frame(Dev,RMSE,RSR,pF,Acc,PCC)
  }
} else {
  rank <- colMeans(rank)
  if (RMSEttest){
    Performance <- data.frame(colMeans(Dev),colMeans(RMSE),colMeans(RSR),colMeans(pF),colMeans(Acc),colMeans(PCC),colMeans(RMSEt))
  } else {
    Performance <- data.frame(colMeans(Dev),colMeans(RMSE),colMeans(RSR),colMeans(pF),colMeans(Acc),colMeans(PCC))
  }
}
Performance <- Performance[rank,]
rownames(Performance) <- Imputations$method
if (RMSEttest) {
  colnames(Performance) <- c('Deviation','RMSE','RSR','p-Value_F-test','Accuracy','PCC','RMSEttest')
} else {
  colnames(Performance) <- c('Deviation','RMSE','RSR','p-Value_F-test','Accuracy','PCC')
}
return(Performance)
}
