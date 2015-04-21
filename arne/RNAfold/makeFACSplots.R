
library(vioplot)


araH <- read.csv(file="araH-FACS.txt",header=FALSE,sep=" ");

araH.DGsort=araH[with(araH,order(V4)), ]
araH.DGsort$V8 = filter(araH.DGsort$V3,rep(1/101,101),sides=2)
pdf("araH-DG.pdf")
plot(araH.DGsort$V4,araH.DGsort$V3,xlab="deltaG",ylab="FACS",main="araH",col="black",type="p",pch = ".")
lines(araH.DGsort$V4,araH.DGsort$V8,col="red")
dev.off()

araH.GCsort=araH[with(araH,order(V5)), ]
araH.GCsort$V8 = filter(araH.GCsort$V3,rep(1/101,101),sides=2)
#pdf("araH-GC.pdf")
#plot(araH.GCsort$V5,araH.GCsort$V3,xlab="%GC",ylab="FACS",main="araH")
#ines(araH.GCsort$V5,araH.GCsort$V8,col="red")
#dev.off()

araH.0<-araH.DGsort$V3[araH.DGsort$V5==0.066667]
araH.1<-araH.DGsort$V3[araH.DGsort$V5==0.066667]
araH.2<-araH.DGsort$V3[araH.DGsort$V5==0.133333]
araH.3<-araH.DGsort$V3[araH.DGsort$V5==0.200000]
araH.4<-araH.DGsort$V3[araH.DGsort$V5==0.266667]
araH.5<-araH.DGsort$V3[araH.DGsort$V5==0.333333]
araH.6<-araH.DGsort$V3[araH.DGsort$V5==0.400000]
araH.7<-araH.DGsort$V3[araH.DGsort$V5==0.466667]
araH.8<-araH.DGsort$V3[araH.DGsort$V5==0.533333]
araH.9<-araH.DGsort$V3[araH.DGsort$V5==0.600000]
araH.10<-araH.DGsort$V3[araH.DGsort$V5==0.666667]
araH.11<-araH.DGsort$V3[araH.DGsort$V5==0.733333]
araH.12<-araH.DGsort$V3[araH.DGsort$V5==0.300000]
araH.13<-araH.DGsort$V3[araH.DGsort$V5==0.366667]
araH.14<-araH.DGsort$V3[araH.DGsort$V5==0.933333]
araH.15<-araH.DGsort$V3[araH.DGsort$V5==1.000000]

#for (i in 1:12 ){
#    araH.[i]<-araH.DGsort$V3[araH.DGsort$V5==i/12]
#}

pdf("araH-GC.pdf")
#plot(araH.GCsort$V5,araH.GCsort$V3,araH.lab="%GC",ylab="FACS",main="narK")
#lines(araH.GCsort$V5,araH.GCsort$V8,col="red")

#vioplot(  araH.3, araH.4 , araH.5 ,araH.6, araH.7,araH.8,araH.9,araH.10,
#        names=c("20.0%","26.7%","33.3%","40.0%","46.7%","53.3%","60.0%","66.7%"),
#        col="grey")
#vioplot(   araH.4 , araH.5 ,araH.6, araH.7,araH.8,araH.9,araH.10,
#        names=c("26.7%","33.3%","40.0%","46.7%","53.3%","60.0%","66.7%"),
#        col="grey")
boxplot(   araH.4 , araH.5 ,araH.6, araH.7,araH.8,araH.9,araH.10,
        names=c("26.7%","33.3%","40.0%","46.7%","53.3%","60.0%","66.7%"),
        col="grey")
title("araH %GC")
dev.off()


pdf("araH-Codon.pdf")
araH.codon1.12<-araH.DGsort$V3[araH.DGsort$V6==12]
araH.codon1.44<-araH.DGsort$V3[araH.DGsort$V6==44]

araH.codon2.33<-araH.DGsort$V3[araH.DGsort$V7==33]
araH.codon2.34<-araH.DGsort$V3[araH.DGsort$V7==34]
araH.codon2.35<-araH.DGsort$V3[araH.DGsort$V7==35]
araH.codon2.36<-araH.DGsort$V3[araH.DGsort$V7==36]
araH.codon2.57<-araH.DGsort$V3[araH.DGsort$V7==57]
araH.codon2.58<-araH.DGsort$V3[araH.DGsort$V7==58]
#vioplot(   araH.codon1.12,araH.codon1.44,
#        araH.codon2.33,araH.codon2.34,araH.codon2.35,araH.codon2.36,araH.codon2.57,araH.codon2.58,
#        names=c("1:ATG","1:ACG","2:TCT","2:TCC","2:TCA","2:TCG","2:AGT","2:AGC"),
#        col="grey")
boxplot(   araH.codon1.12,araH.codon1.44,
        araH.codon2.33,araH.codon2.34,araH.codon2.35,araH.codon2.36,araH.codon2.57,araH.codon2.58,
        names=c("1:ATG","1:ACG","2:TCT","2:TCC","2:TCA","2:TCG","2:AGT","2:AGC"),
        col="grey")

dev.off()

#-----------------------------------------------------------------------
narK <- read.csv(file="narK-FACS.txt",header=FALSE,sep=" ");

narK.DGsort=narK[with(narK,order(V4)), ]
narK.DGsort$V8 = filter(narK.DGsort$V3,rep(1/151,151),sides=2)
pdf("narK-DG.pdf")
plot(narK.DGsort$V4,narK.DGsort$V3,xlab="deltaG",ylab="FACS",main="narK",col="black",type="p",pch = ".")
lines(narK.DGsort$V4,narK.DGsort$V8,col="red")
dev.off()

narK.GCsort=narK[with(narK,order(V5)), ]
narK.GCsort$V8 = filter(narK.GCsort$V3,rep(1/201,201),sides=2)

narK.0<-narK.DGsort$V3[narK.DGsort$V5==0.066667]
narK.1<-narK.DGsort$V3[narK.DGsort$V5==0.066667]
narK.2<-narK.DGsort$V3[narK.DGsort$V5==0.133333]
narK.3<-narK.DGsort$V3[narK.DGsort$V5==0.200000]
narK.4<-narK.DGsort$V3[narK.DGsort$V5==0.266667]
narK.5<-narK.DGsort$V3[narK.DGsort$V5==0.333333]
narK.6<-narK.DGsort$V3[narK.DGsort$V5==0.400000]
narK.7<-narK.DGsort$V3[narK.DGsort$V5==0.466667]
narK.8<-narK.DGsort$V3[narK.DGsort$V5==0.533333]
narK.9<-narK.DGsort$V3[narK.DGsort$V5==0.600000]
narK.10<-narK.DGsort$V3[narK.DGsort$V5==0.666667]
narK.11<-narK.DGsort$V3[narK.DGsort$V5==0.733333]
narK.12<-narK.DGsort$V3[narK.DGsort$V5==0.300000]
narK.13<-narK.DGsort$V3[narK.DGsort$V5==0.366667]
narK.14<-narK.DGsort$V3[narK.DGsort$V5==0.933333]
narK.15<-narK.DGsort$V3[narK.DGsort$V5==1.000000]

#for (i in 1:12 ){
#    narK.[i]<-narK.DGsort$V3[narK.DGsort$V5==i/12]
#}

pdf("narK-GC.pdf")
#plot(narK.GCsort$V5,narK.GCsort$V3,narK.lab="%GC",ylab="FACS",main="narK")
#lines(narK.GCsort$V5,narK.GCsort$V8,col="red")

#vioplot(  narK.3, narK.4 , narK.5 ,narK.6, narK.7,narK.8,narK.9,narK.10,narK.11,
#        names=c("20.0%","26.7%","33.3%","40.0%","46.7%","53.3%","60.0%","66.7%","73.3%"),
#        col="grey")
#vioplot(   narK.4 , narK.5 ,narK.6, narK.7,narK.8,narK.9,narK.10,
#        names=c("26.7%","33.3%","40.0%","46.7%","53.3%","60.0%","66.7%"),
#        col="grey")
boxplot(   narK.4 , narK.5 ,narK.6, narK.7,narK.8,narK.9,narK.10,
        names=c("26.7%","33.3%","40.0%","46.7%","53.3%","60.0%","66.7%"),
        col="grey")
title("narK %GC")

dev.off()
pdf("narK-Codon.pdf")

narK.codon1.33<-narK.DGsort$V3[narK.DGsort$V6==33]
narK.codon1.34<-narK.DGsort$V3[narK.DGsort$V6==34]
narK.codon1.35<-narK.DGsort$V3[narK.DGsort$V6==35]
narK.codon1.36<-narK.DGsort$V3[narK.DGsort$V6==36]
narK.codon1.57<-narK.DGsort$V3[narK.DGsort$V6==57]
narK.codon1.58<-narK.DGsort$V3[narK.DGsort$V6==58]

narK.codon2.21<-narK.DGsort$V3[narK.DGsort$V7==21]
narK.codon2.22<-narK.DGsort$V3[narK.DGsort$V7==22]


#vioplot( narK.codon1.33,narK.codon1.34,narK.codon1.35,narK.codon1.36,narK.codon1.57,narK.codon1.58,
#        narK.codon2.21,narK.codon2.22,
#        names=c("1:TCT","1:TCC","1:TCA","1:TCG","1:AGT","1:AGC","2:CAT","2:CAC"),
#        col="grey")
boxplot( narK.codon1.33,narK.codon1.34,narK.codon1.35,narK.codon1.36,narK.codon1.57,narK.codon1.58,
        narK.codon2.21,narK.codon2.22,
        names=c("1:TCT","1:TCC","1:TCA","1:TCG","1:AGT","1:AGC","2:CAT","2:CAC"),
        col="grey")

dev.off()


data<-read.csv(file="GFP.dat",header=FALSE,sep="\t");
data.DGsort=data[with(data,order(V3)), ]
data.DGsort$V5 = filter(data.DGsort$V4,rep(1/51,51),sides=2)
pdf("gfp.pdf")
plot(data$V3,data$V4,xlab="deltaG",ylab="GFP",main="GFP flourescence",,col="black",type="p",pch = ".")
lines(data.DGsort$V3,data.DGsort$V5,col="red")
dev.off()
