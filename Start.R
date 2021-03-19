Sys.setenv(LANG="en")
library("stats")
library("matrixStats")
source("readInMaxQuant.R")
source("LearnPattern.R")
source("ConstructReferenceData.R")
source("AssignPattern.R")
source("DoImputations.R")

filename <- "TestData.txt"
mtx <- maxQuantToMatrix(filename)
coef <- LearnPattern(mtx)
ref <- ConstructReferenceData(mtx)
sim <- AssignPattern(ref,coef,mtx)

Imputations <- DoImputations(sim)
Performance <- EvaluatePerformance(Imputations,ref,sim)
