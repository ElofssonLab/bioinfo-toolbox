library(vioplot)
genomes=NULL
genomes[1]="escherichia_coli"
genomes[2]="saccharomyces_cerevisae"
genomes[3]="homo_sapiens"
#genome="escherichia_coli"
for (genome in genomes){

file=paste(genome,".df.tsv",sep="")

#dat<-read.table("escherichia_coli.df.tsv", sep='\t', header=T)
#dat<-read.table("saccharomyces_cerevisae.df.tsv", sep='\t', header=T)
#dat<-read.table("homo_sapiens.df.tsv", sep='\t', header=T)
dat<-read.table(file, sep='\t', header=T)



#attach(dat)

dat$Pfam_Meff[dat$Pfam_Meff<1] <- 1
dat$Pfam_Meff[is.na(dat$Pfam_Meff)] <- 1

# fraction of residues with Meff>100 among those that have a valid
# Meff value
length(which(dat$Pfam_Meff > 100))/length(which(dat$Pfam_Meff != -1))
length(which(dat$Pfam_Meff > 100))/length(dat$Pfam_Meff)


# Pfam coverage
length(which(dat$Pfam_ID == ""))/length(dat$aa)

outfile=paste(genome,"-Disorder.pdf",sep="")
pdf(outfile)
hist(dat$Disorder,xlab="Disorder",ylab="Frequency",main=genome)
dev.off()
                                        #hist(dat$Pfam_Meff)
outfile=paste(genome,"-Meff.pdf",sep="")
pdf(outfile)
hist(log10(dat$Pfam_Meff),xlab="log10(Meff)",ylab="Frequency",main=genome)
dev.off()

                                        #plot(Disorder,dat$Pfam_Meff)
#Pfam=dat$Pfam_Meff(replaace

#Fraction disorder with PFam hit
loop=seq(1,10)
disoPfam=NULL
disofrac=NULL
Meff100=NULL
Meff200=NULL
Meff500=NULL
Meff1000=NULL
frac=NULL
for (i in loop){
    frac[i]=0.1*i
    disofrac[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))/length(dat$Disorder)
    disoPfam[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1  & dat$Pfam_ID != ""))/length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))
    Meff100[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1  & dat$Pfam_Meff >100 ))/length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))
    Meff200[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1  & dat$Pfam_Meff >200 ))/length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))
    Meff500[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1  & dat$Pfam_Meff >500 ))/length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))
    Meff1000[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1  & dat$Pfam_Meff >1000 ))/length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))
}


ntrees <- 10
linetype <- c(1:ntrees)



outfile=paste(genome,"-Meff.pdf",sep="")
pdf(outfile)
plot(frac,disoPfam,type="l",ylim=c(0,1),ylab="%residues",xlab="% disorder",main=genome,lty=linetype[1])
points(frac,disofrac,col="red",lty=linetype[1],type="p")
lines(frac,Meff100,col="grey",lty=linetype[2])
lines(frac,Meff200,col="grey",lty=linetype[3])
lines(frac,Meff500,col="grey",lty=linetype[4])
lines(frac,Meff1000,col="grey",lty=linetype[5])
legend(0.8,0.8, c("Pfam","100","200","500","1000"),lty=linetype,cex=0.5)
dev.off()

    
disoMeff0<-dat$Pfam_Meff[dat$Disorder > 0][dat$Disorder <0.1]
disoMeff1<-dat$Pfam_Meff[dat$Disorder > 0.1][dat$Disorder <0.2]
disoMeff2<-dat$Pfam_Meff[dat$Disorder > 0.2][dat$Disorder <0.3]
disoMeff3<-dat$Pfam_Meff[dat$Disorder > 0.3][dat$Disorder <0.4]
disoMeff4<-dat$Pfam_Meff[dat$Disorder > 0.4][dat$Disorder <0.5]
disoMeff5<-dat$Pfam_Meff[dat$Disorder > 0.5][dat$Disorder <0.6]
disoMeff6<-dat$Pfam_Meff[dat$Disorder > 0.6][dat$Disorder <0.7]
disoMeff7<-dat$Pfam_Meff[dat$Disorder > 0.7][dat$Disorder <0.8]
disoMeff8<-dat$Pfam_Meff[dat$Disorder > 0.8][dat$Disorder <0.9]
disoMeff9<-dat$Pfam_Meff[dat$Disorder > 0.9][dat$Disorder <1.0]

outfile=paste(genome,"-DisoBox.pdf",sep="")
pdf(outfile)
boxplot(disoMeff0,disoMeff1,disoMeff2,disoMeff3,disoMeff4,disoMeff5,disoMeff6,disoMeff7,disoMeff8,disoMeff9,names=c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95),log="y",main=genome,ylab="%residues",xlab="% disorder")
dev.off()


                                        # TM (not important)

TMvar=NULL
TMvar[1]=1
TMvar[2]=2
TMvar[3]=3

TMPFam=NULL
TMPFam[1]=length(which(dat$TM_topogy== "M"  & dat$Pfam_ID != "" ))/length(which(dat$TM_topogy== "M"))
TMPFam[2]=length(which(dat$TM_topogy== "i"  & dat$Pfam_ID != "" ))/length(which(dat$TM_topogy== "i"))
TMPFam[3]=length(which(dat$TM_topogy== "o"  & dat$Pfam_ID != "" ))/length(which(dat$TM_topogy== "o"))

TMMeff100=NULL
TMMeff100[1]=length(which(dat$TM_topogy== "M"  & dat$Pfam_Meff >100 ))/length(which(dat$TM_topogy== "M"))
TMMeff100[2]=length(which(dat$TM_topogy== "i"  & dat$Pfam_Meff >100 ))/length(which(dat$TM_topogy== "i"))
TMMeff100[3]=length(which(dat$TM_topogy== "o"  & dat$Pfam_Meff >100 ))/length(which(dat$TM_topogy== "o"))

TMMeff200=NULL
TMMeff200[1]=length(which(dat$TM_topogy== "M"  & dat$Pfam_Meff >200 ))/length(which(dat$TM_topogy== "M"))
TMMeff200[2]=length(which(dat$TM_topogy== "i"  & dat$Pfam_Meff >200 ))/length(which(dat$TM_topogy== "i"))
TMMeff200[3]=length(which(dat$TM_topogy== "o"  & dat$Pfam_Meff >200 ))/length(which(dat$TM_topogy== "o"))

TMMeff500=NULL
TMMeff500[1]=length(which(dat$TM_topogy== "M"  & dat$Pfam_Meff >500 ))/length(which(dat$TM_topogy== "M"))
TMMeff500[2]=length(which(dat$TM_topogy== "i"  & dat$Pfam_Meff >500 ))/length(which(dat$TM_topogy== "i"))
TMMeff500[3]=length(which(dat$TM_topogy== "o"  & dat$Pfam_Meff >500 ))/length(which(dat$TM_topogy== "o"))

TMMeff1000=NULL
TMMeff1000[1]=length(which(dat$TM_topogy== "M"  & dat$Pfam_Meff >1000 ))/length(which(dat$TM_topogy== "M"))
TMMeff1000[2]=length(which(dat$TM_topogy== "i"  & dat$Pfam_Meff >1000 ))/length(which(dat$TM_topogy== "i"))
TMMeff1000[3]=length(which(dat$TM_topogy== "o"  & dat$Pfam_Meff >1000 ))/length(which(dat$TM_topogy== "o"))

outfile=paste(genome,"-TM.pdf",sep="")
pdf(outfile)
plot(TMvar,TMPFam,type="p",ylim=c(0,1),ylab="%residues",xlab="TM",names=c("M","i","o"))
lines(TMvar,TMMeff100,col="blue")
lines(TMvar,TMMeff200,col="red")
lines(TMvar,TMMeff500,col="green")
lines(TMvar,TMMeff1000,col="grey")
dev.off()
}    
#Pfam.log <- log10(dat$Pfam_Meff)
