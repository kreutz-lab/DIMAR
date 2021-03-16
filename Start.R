Sys.setenv(LANG="en")
library("stats")
source("readInMaxQuant.R")
source("LearnPattern.R")
source("ConstructReferenceData.R")
source("AssignPattern.R")

filename <- "TestData.txt"
mtx <- maxQuantToMatrix(filename)
mtx <- mtx[1:100,1:10]

fit <- LearnPattern(mtx)
sim <- ConstructReferenceData(mtx)
sim <- AssignPattern(sim,fit)
