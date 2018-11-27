#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
infile=args[1]
outfile=args[2]
d<-read.csv(file=infile,sep=" ",header=F)
d$V5<-scale(d$V4)
write.table(d,file=outfile,sep=" ",row.names=FALSE,col.names=FALSE)
