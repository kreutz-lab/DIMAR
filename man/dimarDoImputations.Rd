% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimarDoImputations.R
\name{dimarDoImputations}
\alias{dimarDoImputations}
\title{dimarDoImputations}
\usage{
dimarDoImputations(mtx, method = "fast", lib = NULL)
}
\arguments{
\item{mtx}{Quantitative matrix}

\item{method}{Imputation method(s)}

\item{lib}{R packages of imputation methods (to be loaded in parallel loop 'foreach')}
}
\value{
Imputed data set(s)
}
\description{
Performs imputations
}
\examples{
mtx <- matrix(rnorm(1000),nrow=100)
mtx[sample(c(1:1000),100)] <- NA
Imputations <- dimarDoImputations(mtx, c('impSeqRob', 'ppca', 'imputePCA'))
}
