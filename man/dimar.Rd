% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimar.R
\name{dimar}
\alias{dimar}
\title{dimar}
\usage{
dimar(mtx, pattern = NULL, methods = "fast", npat = NULL, group = "cluster")
}
\arguments{
\item{mtx}{Data matrix or MaxQuant input file ('.txt')}

\item{pattern}{Search pattern for specifying sample names read in as default data, if not specified the user will be asked}

\item{methods}{List of imputation algorithms ['fast'] uses the nine most selected algorithms}

\item{npat}{Number of patterns of MVs to be simulated and to test the algorithms on [5/10/20 depending on the size of the data]}

\item{group}{vector of group indices for ttest (group==1 vs group==2) ['cluster'] as a default, clustering with 2 cluster is performed}
}
\description{
DIMAR: Data-driven selection of an imputation algorithm in R
}
\examples{
mtx <- matrix(rnorm(1000), nrow=100)
mtx[sample(c(1:1000),100)] <- NA
Imp <- dimar(mtx)

filename <- "Test1.txt"
filepath <- system.file("extdata", filename, package = "DIMAR")
Imp <- dimar(mtx = filepath, pattern = 'Intensity', group = c('PKB','PKC'))
}
