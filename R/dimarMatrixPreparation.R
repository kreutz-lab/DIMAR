#' dimarMatrixPreparation
#' 
#' @description Prepare/Filter matrix
#' @param mtx Quantitative matrix
#' @param nacut minimum number of measured data points
#' @param logflag should logarithm be applied (true,false,log2,log10)
#' @export dimarMatrixPreparation dimarMatrixPreparation
#' @examples
#' mtx <- dimarMatrixPreparation(mtx)

dimarMatrixPreparation <- function(mtx,nacut=1,logflag='auto'){

    # if no NA: 0 -> NA
    mtx[mtx == "NaN"] <- NA
    if (sum(is.na(mtx))==0) {
        mtx[mtx == 0] <- NA
    }
    if (sum(is.na(mtx))==0) {
        print('No MVs in your dataset.')
    }
    
    # logarithm
    switch(logflag,
           auto={ if (max(mtx,na.rm=T)>100) {
                    mtx <- log2(mtx) } },
           true={mtx <- log2(mtx)},
           log2={mtx <- log2(mtx)},
           log10={mtx <- log10(mtx)},
           warning(paste('dimarMatrixPreparation: logflag',logflag,'not known. Expand code here. No transformation is performed.'))
    )
    
    # Nacut
    if (nacut>=0 && nacut<1) {
        mtx = mtx[rowSums(!is.na(mtx))>nacut*dim(mtx)[2],]
    } else if (nacut>=1){
        mtx = mtx[rowSums(!is.na(mtx))>nacut,]
    } else {
        warning(paste('dimarMatrixPreparation: nacut',nacut,'not known. Expand code here. No transformation is performed.'))
    }
    
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