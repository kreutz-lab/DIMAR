% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimarEvaluatePerformance.R
\name{dimarEvaluatePerformance}
\alias{dimarEvaluatePerformance}
\title{dimarEvaluatePerformance}
\usage{
dimarEvaluatePerformance(
  Imputations,
  ref,
  sim,
  rankby = "RMSE",
  RMSEttest = TRUE,
  group = "cluster"
)
}
\arguments{
\item{Imputations}{Imputed data set(s)}

\item{ref}{Reference data}

\item{sim}{Simulated patterns of MVs}

\item{rankby}{Performance measure which should serve as rank criterion}

\item{RMSEttest}{flag if RMSE of ttest should be calculated}

\item{group}{indices for ttest}
}
\value{
Data frame containing the following performance measures for
each imputation method: Deviation, RMSE, RSR, p-Value_F-test, Accuracy, PCC, and in case
of RMSEttest=TRUE the RMSE t-test result
}
\description{
Evaluates performance of imputation algorithms.
}
\examples{
mtx <- matrix(rnorm(1000),nrow=100)
mtx[sample(c(1:1000),100)] <- NA
coef <- dimarLearnPattern(mtx)
ref <- dimarConstructReferenceData(mtx)
sim <- dimarAssignPattern(ref, coef, mtx)
Imputations <- dimarDoImputations(sim, c('impSeqRob', 'ppca', 'imputePCA'))
Performance <- dimarEvaluatePerformance(Imputations, ref, sim)
}
