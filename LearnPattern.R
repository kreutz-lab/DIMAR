LearnPattern <- function(mtx) {
  
# Subsample indices
nsub <- ceiling(dim(mtx)[1]/5)
npersub <- ceiling(dim(mtx)[1]/nsub)
indrand <- sample(1:dim(mtx)[1],dim(mtx)[1])

for (i in 1:nsub) {
  #  i <- 1
  ind <- indrand[(npersub*(i-1)+1):(npersub*i)]
  design <- ConstructDesignMatrix(mtx[ind,])
  design <- ConstructRegularizationMatrix(design)
  
  #fit <- glm.fit(X,y,family=binomial(),weights=rep(1,dim(X)[1]))
  fit <- glm.fit(design$X,design$y, family=binomial())
  if (i==1) {
    coef <- coefficients(fit)
  } else {
    coef <- rbind(coef,coefficients(fit))
  }
}
# row coefficients stay same, intensity/column coefficients are set to mean over nsub (for loop)
coefsubmean <- c(colMeans(coef[,design$Xtype!=3]),coef[,design$Xtype==3])
return(coefsubmean)
}