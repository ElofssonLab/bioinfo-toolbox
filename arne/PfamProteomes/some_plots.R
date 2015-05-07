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


# Pie Charts

# Fraction of residues with Pfam-coverage and Meff coverage=NULL
coverage=NULL
labels=NULL
coverage[1]=length(which(dat$Pfam_ID == ""))
labels[1]="NoPfam"
coverage[2]=length(which(dat$Pfam_Meff > 1000))
labels[2]=">1000"
coverage[3]=length(which(dat$Pfam_Meff > 500 & dat$Pfam_Meff <= 1000))
labels[3]=">500"
coverage[4]=length(which(dat$Pfam_Meff > 200 & dat$Pfam_Meff <= 500))
labels[4]=">200"
coverage[5]=length(which(dat$Pfam_Meff > 100 & dat$Pfam_Meff <= 200))
labels[5]=">100"
coverage[6]=length(which(dat$Pfam_ID != "" & dat$Pfam_Meff))
labels[6]="Pfam"
pct <- round(coverage/sum(coverage)*100)
outfile=paste(genome,"-coverage.pdf",sep="")
pdf(outfile)
pie(coverage,labels = labels, col=rainbow(length(labels)),main=genome)
dev.off()

# FOr ordered and disordered residues
covorder=NULL
labels=NULL
colors=NULL
covorder[1]=length(which(dat$Pfam_ID == ""))
labels[1]="NoPfam"
colors[1]="cadetblue1"
covorder[6]=length(which(dat$Pfam_Meff > 1000 & dat$Disorder<0.5))
labels[6]=">1000"
colors[6]="blue4"
covorder[5]=length(which(dat$Pfam_Meff > 500 & dat$Pfam_Meff <= 1000 & dat$Disorder<0.5))
labels[5]=">500"
colors[5]="blue3"
covorder[4]=length(which(dat$Pfam_Meff > 200 & dat$Pfam_Meff <= 500 & dat$Disorder<0.5))
labels[4]=">200"
colors[4]="blue2"
covorder[3]=length(which(dat$Pfam_Meff > 100 & dat$Pfam_Meff <= 200 & dat$Disorder<0.5))
labels[3]=">100"
colors[3]="blue1"
covorder[2]=length(which(dat$Pfam_ID != "" & dat$Pfam_Meff & dat$Disorder<0.5))
labels[2]="Pfam"
colors[2]="cadetblue2"
pct <- round(covorder/sum(covorder)*100)
labels <- paste(labels,pct)
labels <- paste(labels,"%",sep="")
outfile=paste(genome,"-pie-order.pdf",sep="")
pdf(outfile)
pie(covorder,labels = labels, col=rainbow(length(labels)),main=genome)
dev.off()

covdisorder=NULL
dislabels=NULL
covdisorder[1]=length(which(dat$Pfam_ID == ""))
covorder[7]=length(which(dat$Pfam_ID == ""))
dislabels[1]="NoPfam"
labels[7]="Diso-NoPfam"
colors[7]="lightpink"

covdisorder[6]=length(which(dat$Pfam_Meff > 1000 & dat$Disorder>=0.5))
covorder[12]=length(which(dat$Pfam_Meff > 1000 & dat$Disorder>=0.5))
dislabels[6]=">1000"
labels[12]=">1000 diso"
colors[12]="red"
covdisorder[5]=length(which(dat$Pfam_Meff > 500 & dat$Pfam_Meff <= 1000 & dat$Disorder>=0.5))
covorder[11]=length(which(dat$Pfam_Meff > 500 & dat$Pfam_Meff <= 1000 & dat$Disorder>=0.5))
dislabels[5]=">500"
labels[11]=">500 diso"
colors[9]="deepink3"
covdisorder[4]=length(which(dat$Pfam_Meff > 200 & dat$Pfam_Meff <= 500 & dat$Disorder>=0.5))
covorder[10]=length(which(dat$Pfam_Meff > 200 & dat$Pfam_Meff <= 500 & dat$Disorder>=0.5))
dislabels[4]=">200"
labels[10]=">200 diso"
colors[9]="deepink2"
covdisorder[3]=length(which(dat$Pfam_Meff > 100 & dat$Pfam_Meff <= 200 & dat$Disorder>=0.5))
covorder[9]=length(which(dat$Pfam_Meff > 100 & dat$Pfam_Meff <= 200 & dat$Disorder>=0.5))
dislabels[3]=">100"
labels[9]=">100 diso"
colors[9]="deeppink1"
covdisorder[2]=length(which(dat$Pfam_ID != "" & dat$Pfam_Meff & dat$Disorder>=0.5))
covorder[8]=length(which(dat$Pfam_ID != "" & dat$Pfam_Meff & dat$Disorder>=0.5))
dislabels[2]="Pfam"
labels[8]="Pfam diso"
colors[8]="lightpink2"

pct <- round(covdisorder/sum(covdisorder)*100,digits=2)
dislabels <- paste(dislabels,pct)
dislabels <- paste(dislabels,"%",sep="")
outfile=paste(genome,"-pie-disorder.pdf",sep="")
pdf(outfile)
pie(covdisorder,labels = labels, col=rainbow(length(labels)),main=genome)
dev.off()

pct <- round(covorder/sum(covorder)*100,digits=2)
labels <- paste(labels,pct)
labels <- paste(labels,"%",sep="")
outfile=paste(genome,"-pie-all.pdf",sep="")
pdf(outfile)
pie(covorder,labels = labels, col=colors   ,main=genome,cex=0.3)
dev.off()



}    
#Pfam.log <- log10(dat$Pfam_Meff)
