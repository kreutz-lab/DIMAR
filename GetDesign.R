GetDesign <- function(mtx) {

X <- matrix(0L,nrow=dim(mtx)[1]*dim(mtx)[2],ncol=dim(mtx)[1]+dim(mtx)[2]+1)
X[,1] <- scale(rep(rowMeans(mtx,na.rm=TRUE),dim(mtx)[2]))
Xname <- 'mean'
Xtype <- 1
row <- rep(1:dim(mtx)[1],dim(mtx)[2])
col <- rep(1:dim(mtx)[2],each=dim(mtx)[1])
for (i in 1:dim(mtx)[1]) {
  X[col==i,i+1] <- 1
  Xname <- c(Xname,paste('Col',i))
  Xtype <- c(Xtype,2)
}
for (i in 1:dim(mtx)[2]) {
  X[row==i,i+1+dim(mtx)[1]] <- 1
  Xname <- c(Xname,paste('Row',i))
  Xtype <- c(Xtype,3)
}
y <- as.numeric(is.na(mtx))
#return(list(y,X,Xtype))
return(list(X,Xtype))
}