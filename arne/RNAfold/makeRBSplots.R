library(vioplot)

narK<-read.csv(file="narK-RBS.out",header=FALSE,sep=" ");
araH<-read.csv(file="araH-RBS.out",header=FALSE,sep=" ");

araH$V5 = log10(araH$V3)
narK$V5 = log10(narK$V3)

araH.sort=araH[with(araH,order(V5)), ]
narK.sort=narK[with(narK,order(V5)), ]

araH.sort$V6 = filter(araH.sort$V4,rep(1/101,101),sides=2)
narK.sort$V6 = filter(narK.sort$V4,rep(1/101,101),sides=2)




pdf("araH-RBS.pdf")
plot (araH.sort$V5,araH.sort$V4,col="grey",xlab="log10(TIR)",ylab="FACS",type="p",cex=0.5)
lines(araH.sort$V5,araH.sort$V6,col="red",lwd=4)
dev.off()
pdf("narK-RBS.pdf")
plot (narK.sort$V5,narK.sort$V4,col="grey",xlab="log10(TIR)",ylab="FACS",type="p",cex=0.5)
lines(narK.sort$V5,narK.sort$V6,col="red",lwd=4)
dev.off()

pdf("araH-RBS-box.pdf")
araH.0<-araH$V3[araH$V4>=4.0 ][ araH$V4<4.2 ]
araH.1<-araH$V3[araH$V4>=4.2 ][ araH$V4<4.4 ]
araH.2<-araH$V3[araH$V4>=4.4 ][ araH$V4<4.6 ]
araH.3<-araH$V3[araH$V4>=4.6 ][ araH$V4<4.8 ]
araH.4<-araH$V3[araH$V4>=4.8 ][ araH$V4<5.0 ]
boxplot(araH.0,araH.1,araH.2,araH.3,araH.4,names=c(4.1,4.3,4.5,4.7,4.9),log="y")
dev.off()

pdf("narK-RBS-box.pdf")
narK.0<-narK$V3[narK$V4>=4.0 ][ narK$V4<4.2 ]
narK.1<-narK$V3[narK$V4>=4.2 ][ narK$V4<4.4 ]
narK.2<-narK$V3[narK$V4>=4.4 ][ narK$V4<4.6 ]
narK.3<-narK$V3[narK$V4>=4.6 ][ narK$V4<4.8 ]
narK.4<-narK$V3[narK$V4>=4.8 ][ narK$V4<5.0 ]
boxplot(narK.0,narK.1,narK.2,narK.3,narK.4,names=c(4.1,4.3,4.5,4.7,4.9),log="y")
dev.off()

