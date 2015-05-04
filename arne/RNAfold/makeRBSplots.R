library(vioplot)

narK<-read.csv(file="narK-RBS.out",header=FALSE,sep=" ");
araH<-read.csv(file="araH-RBS.out",header=FALSE,sep=" ");

araH.sort=araH[with(araH,order(V3)), ]
narK.sort=narK[with(narK,order(V3)), ]



pdf("araH-RBS.pdf")
plot(araH$V4,araH$V3,col="black",ylab="TIR",xlab="FACS",log="y")
araH.fit =lm (araH$V3~log(araH$V4))
abline(araH.fit,col="blue",lwd=4)
dev.off()

pdf("narKH-RBS.pdf")
 plot(narK$V4,narK$V3,col="black",ylab="TIR",xlab="FACS",log="y")
dev.off()

pdf("araH-RBS2.pdf")
plot (narK$V3,narK$V4,col="red",ylab="TIR",xlab="FACS",log="x")
dev.off()
pdf("araH-RBS2.pdf")
plot (narK$V3,narK$V4,col="red",ylab="TIR",xlab="FACS",log="x")
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

