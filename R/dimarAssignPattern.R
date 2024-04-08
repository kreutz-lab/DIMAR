#' dimarAssignPattern
#'
#' @description Assigns patterns of MVs to reference data set with a logistic regression model
#' @return Reference dataset with applied pattern
#' @param ref Reference dataset
#' @param coef Logistic regression coefficients of missing value pattern
#' @param mtx Quantitative matrix
#' @param npat Number of patterns
#' @param dontAllowOnlyNA Default: TRUE ; Proteins with only NA will be assigned a random value. 
#' If DIMAR is used to simulate by group, set this to FALSE. Proteins with only missing values in one of the groups are common and important for a realistic missing value pattern.
#' @export dimarAssignPattern
#' @examples
#' mtx <- matrix(rnorm(1000),nrow=100)
#' mtx[sample(c(1:1000),100)] <- NA
#' coef <- dimarLearnPattern(mtx)
#' ref <- dimarConstructReferenceData(mtx)
#' sim <- dimarAssignPattern(ref, coef, mtx)

dimarAssignPattern <- function(ref, coef, mtx = NULL, npat = NULL, dontAllowOnlyNA = TRUE, seed = NULL, groupDesign = rep(c(1,2), each = ncol(ref)/2)) {
  
  #if a seed is given, draw from coefs the same way as coefficients from potential previous linear model intensity simulation (lmInt-Sim) were drawn. 
  #In order to conserve missing value pattern it is necessary to ensure features are simulated with intensity and missing value coefficients estimated from the same feature in the input.
  if (!is.null(seed))
    coef <- drawCoefWithSeed(ref = ref, coef = coef, seed = seed, groupDesign = groupDesign)
  
  if (is.null(npat)) {
    if (nrow(ref)*ncol(ref) < 50000) {
      npat <- 20
    } else if (nrow(ref)*ncol(ref) < 100000) {
      npat <- 10
    } else{
      npat <- 5
    }
  }

  X <- dimarConstructDesignMatrix(ref)
  pat <- array(NA,c(dim(ref),npat))

  for (i in 1:npat) {
    # probability of MV by logistic regression
    yhat <- exp(X$X %*% coef)/(1 + exp(X$X %*% coef))
    p <- matrix(yhat,nrow = nrow(ref))
    # binomial draw
    r <- matrix(stats::runif(length(yhat)),nrow = dim(p)[1])
    ind <- which(r < p & !is.na(ref),arr.ind = TRUE)
    # if ref has MV already, #MV = #MV of original mtx
    if (!is.null(mtx) & nrow(ind) > sum(is.na(mtx))) {
      ind <- ind[sample(1:nrow(ind),sum(is.na(mtx))),]
    }
    # Assign NA to ref
    pat1 <- ref
    pat1[ind] <- NA
    if(dontAllowOnlyNA){
    # if protein not measured at all, randomly assign one data point
      if (any(rowSums(is.na(pat1)) == ncol(pat1))) {
        allna <- which(rowSums(is.na(pat1)) == ncol(pat1))
        allna <- cbind(allna,sample(1:ncol(pat1),length(allna),replace=TRUE))
        pat1[allna] <- ref[allna]
      }
    }
    pat[,,i] <- pat1
  }
  print(paste(npat,'patterns of MVs are assigned.'))
  return(pat)
  
  #The Idea for this function is to ensure that the missing values and the intensity value of a certain feature are simulated with a coefficients estimated from the same feature of the input matrix. 
  #It draws the coefficients estimated from the Logistic Regression within DIMAR with the same seed as the coefficients estimated for the Linear Model were previously drawn with.
  drawCoefWithSeed <- function(ref, coef, seed = NULL, groupDesign){  
    DE_idx <- grep("DE",rownames(ref))
    colCoef <- coef[grep("col",names(coef))]
    DECoef <- coef[grep("DE",names(coef))]
    rowCoef <- coef[grep("row", names(coef))]
    rowCoef <- rowCoef[!grepl("DE",names(rowCoef))] #exclude rows that are differentially expressed, as values for DE proteins are drawn separately
    set.seed(seed)
    newColCoef <- sample(colCoef, ncol(ref), replace = T) #draw column coefficients for each column present in ref
    allRowCoef <- c()
    for(g in unique(groupDesign)){
      DECoef_group <- DECoef[grep(paste0("#",g), names(DECoef))]
      rowCoef_group <- rowCoef[grep(paste0("#",g),names(rowCoef))]
      set.seed(seed)
      drawn_rowCoef_group <- sample(rowCoef_group, nrow(ref), replace = T) #first draw for each feature a coefficient from not DE features
      names(drawn_rowCoef_group) <- sub("#[0-9]","",names(drawn_rowCoef_group)) #remove group assignment
      set.seed(seed)
      drawn_rowCoef_group[DE_idx] <- sample(DECoef_group, length(DE_idx), replace = T) #replace coefficients for DE features with coefficients estimated from DE features
      names(drawn_rowCoef_group)[DE_idx] <- paste0(names(drawn_rowCoef_group[DE_idx]), "_DE") #add DE label to name
      names(drawn_rowCoef_group) <- paste0(names(drawn_rowCoef_group),"#",g) #reintroduce group label
      drawn_rowCoef_group
      allRowCoef <- c(allRowCoef, drawn_rowCoef_group)
    }
    newCoef <- c(coef[1:2],newColCoef,allRowCoef)
    
    return(newCoef)
  }
  }  

