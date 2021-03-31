
dimarMatrixPreparation <- function(dat){

    mtx[mtx == "NaN"] <- NA
    if (sum(is.na(mtx))==0) {
        mtx[mtx == 0] <- NA
    }
    if (sum(is.na(mtx))==0) {
        print('No MVs in your dataset.')
    }
    
    # log2 transform intensities
    mtx <- log2(mtx)
    
    # check for proteins with NA intensity across all samples
    allColNA <- as.vector(apply(mtx, 1, function(r) {
        return(all(is.na(r)))
    }))
    mtx <- mtx[!allColNA, ]
    message(paste("Number of proteins with empty entries:",
        length(which(allColNA))))
    
    mtx<- `dimnames<-`(`dim<-`(as.numeric(mtx), dim(mtx)), dimnames(mtx))
    
    return(mtx)
}