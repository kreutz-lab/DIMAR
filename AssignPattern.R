AssignPattern <- function(mtx,fit,npat=0) {

  if (npat==0){
    if (dim(mtx)[1]*dim(mtx)[2]<50000) {
      npat<-20
    } else if (dim(mtx)[1]*dim(mtx)[2]<100000){
      npat<-10
    } else{
      npat<-5
    }
  }
  npat <-1
  X <- ConstructDesignMatrix(mtx)
  sim <- matrix(NA,dim(mtx)[1],dim(mtx)[2])
  for (i in 1:npat){
    yhat <- predict.glm(fit)
  }
  return(yhat)
}