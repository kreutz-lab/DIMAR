#' dimarReadInMaxQuant
#'
#' @description Transforms MaxQuant proteinGroups.txt file content to quantitative matrix
#' @return Quantitative matrix
#' @param filename String of MaxQuant proteinGroups file name
#' @export dimarReadInMaxQuant
#' @examples
dimarReadInMaxQuant <- function(filename){
    # Read file
    dat <- read.csv(filename, allowEscapes = TRUE, check.names = FALSE,sep = "\t")

    # Select all columns with LFQ intensities
    mtx <- as.matrix(dat[, grepl("^Int", names(dat))])
    mtx[mtx == 0] <- NA
    mtx[mtx == "NaN"] <- NA
    # convert character matrix to numeric matrix, strings that are possibly
    # left in matrix by accident are converted to NAs
    mtx <- `dimnames<-`(`dim<-`(as.numeric(mtx), dim(mtx)), dimnames(mtx))

    # check for proteins with NA intensity across all samples
    allColNA <- as.vector(apply(mtx, 1, function(r) {
        return(all(is.na(r)))
    }))
    message(paste("Number of proteins with empty entries:",
                  length(which(allColNA))))

    if (!("Protein IDs" %in% colnames(dat))) {
        colnames(dat)[1] <- "Protein IDs"
    }
    # check if exist and append to featureAnnotations
    featureAnnotations <- data.frame(proteinName = dat[, "Protein IDs"])
    annotations <- c(
        "Fasta headers", "Q-value", "Reverse", "Peptides", "Reverse",
        paste(c("Potential contaminant","Contaminant"),collapse = "|"),
        "Only identified by site")
    fieldnames <- c("proteinDescription", "idScore","isDecoy", "nbPeptides",
                    "isFiltered", "isPotential.contaminant", "isIdentified.by.site")

    # check for potential contaminant and only identfied by side proteins
    bool1 <- dat[,grep(annotations[6],colnames(dat),value = TRUE)] == "+" & !is.na(
        dat[,grep(annotations[6],colnames(dat),value = TRUE)])
    bool2 <- dat[["Only identified by site"]] == "+" & !is.na(
        dat[["Only identified by site"]])
    bool3 <- dat[["Reverse"]] == "+" & !is.na(dat[["Reverse"]])
    # logical is better because it has same dimension as the data
    ixs <-  bool1 | bool2 | bool3

    for (i in seq_len(length(fieldnames))) {
        if (length(grep(annotations[i], colnames(dat))) > 0) {
            if (length(grep(fieldnames[i], c("isDecoy",
                                           "isPotential.contaminant",
                                           "Only identified by site"))) > 0) {
                featureAnnotations[[fieldnames[i]]] <-
                    dat[, grep(annotations[i], colnames(dat), value = TRUE)] == "+"
            } else {
                featureAnnotations[[fieldnames[i]]] <- dat[,annotations[i]]
            }
        }
    }

    # featureAnnotations <- featureAnnotations[!allColNA, ]
    row.names(mtx) <- dat[, "Protein IDs"]

    # # remove empty rows
    mtx <- as.matrix(mtx[!allColNA, ])

    # log2 transform intensities
    mtx <- log2(mtx)

  #  featureAnnotations <- as.list(featureAnnotations)
  #  featureAnnotations[["ixs"]]<-ixs

  #  sexp <- SummarizedExperiment(assays=list(data=mtx),
  #                               rowData=featureAnnotations,
  #                               metadata=list(pxdid = filename))
    # remove empty rows
  #  sexp <- sexp[!allColNA, ]
    # write.csv(mtx[!allColNA, ], file=paste0(filename, ".csv"))
    return(mtx)
}
