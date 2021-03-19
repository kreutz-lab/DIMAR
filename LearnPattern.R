source("ConstructDesignMatrix.R")
source("ConstructRegularizationMatrix.R")

LearnPattern <- function(mtx) {
  
# Subsample indices
if (dim(mtx)[1]>1000) {
  nsub <- ceiling(dim(mtx)[1]/1000)
  indrand <- sample(1:dim(mtx)[1],dim(mtx)[1])
  npersub <- ceiling(dim(mtx)[1]/nsub)
} else { 
  nsub<-1 
  ind <- 1:dim(mtx)[1]
}

for (i in 1:nsub) {
  #  i <- 1
  if (nsub>1){
    ind <- indrand[(npersub*(i-1)+1):(npersub*i)]
  }
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
# sort row coefficients, intensity/column coefficients are set to mean over nsub (for loop)
if (nsub>1){
  coef <- c(colMeans(coef[,design$Xtype!=3]),sort(coef[,design$Xtype==3]))
}

  print('Pattern of MVs is learned by logistic regression.')
  return(coef)
}