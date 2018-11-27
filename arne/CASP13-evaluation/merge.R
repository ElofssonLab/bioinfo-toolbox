#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
infile1=args[1]
infile2=args[2]
infile3=args[3]
outfile=args[4]
a<-read.csv(file=infile1,sep=" ",header=F)
b<-read.csv(file=infile2,sep=" ",header=F)
c<-read.csv(file=infile3,sep=" ",header=F)
temp<-merge(a,b,by="V1", sort = TRUE)
final<-merge(temp,c,by="V1", sort = TRUE)
head(final)
final$average<-(final$V5.x+final$V5.y+final$V5)/3
sorted<-final[order(-final$average),] 
write.table(sorted,file=outfile,sep=" ",row.names=FALSE,col.names=FALSE)
