
library(DIMAR)
filename <- "TestData.txt"
filepath <- system.file("extdata", filename, package = "DIMAR") 
Imp <- dimar(filepath)
