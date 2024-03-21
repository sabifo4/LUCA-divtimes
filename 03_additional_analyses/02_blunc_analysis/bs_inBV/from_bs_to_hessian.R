#! /bin/env Rscript


############################################
library(MASS)


############################################
args = commandArgs(trailingOnly=TRUE)

infile = args[1]
outfile = args[2]


############################################
m <- as.matrix(read.table(infile))

c <- cov(m)

c <- c + diag(ncol(m))*1e-6

h <- -solve(c)

write.table(h, file=outfile, col.names=FALSE, row.names=F)


