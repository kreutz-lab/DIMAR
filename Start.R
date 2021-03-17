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
mtx <- mtx[1:100,1:10]

fit <- LearnPattern(mtx)
ref <- ConstructReferenceData(mtx)
sim <- AssignPattern(ref,fit)

imp <- DoImputations(mtx)
