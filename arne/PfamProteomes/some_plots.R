library(vioplot)

                                        #dat<-read.table("escherichia_coli.df.tsv", sep='\t', header=T)
#dat<-read.table("saccharomyces_cerevisae.df.tsv", sep='\t', header=T)
dat<-read.table("homo_sapiens.df.tsv", sep='\t', header=T)


attach(dat)

Pfam_Meff[Pfam_Meff==-1] <- 1

# fraction of residues with Meff>100 among those that have a valid
# Meff value
length(which(Pfam_Meff > 100))/length(which(Pfam_Meff != -1))
length(which(Pfam_Meff > 100))/length(Pfam <- meff)


# Pfam coverage
length(which(Pfam_ID == ""))/length(aa)

hist(Disorder)
hist(Pfam_Meff)
hist(log10(Pfam_Meff))
#plot(Disorder,Pfam_Meff)
Pfam=Pfam_Meff(replaace

#Fraction disorder with PFam hit
loop=seq(1,10)
disofrac=NULL
Meff100=NULL
Meff200=NULL
Meff500=NULL
    Meff1000=NULL
    frac=NULL
    for (i in loop){
        frac[i]=0.1*i
        disofrac[i]=length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1  & Pfam_ID != ""))/length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1))
        Meff100[i]=length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1  & Pfam_Meff >100 ))/length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1))
        Meff200[i]=length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1  & Pfam_Meff >200 ))/length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1))
        Meff500[i]=length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1  & Pfam_Meff >500 ))/length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1))
        Meff1000[i]=length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1  & Pfam_Meff >1000 ))/length(which(Disorder > (i-1)*0.1 & Disorder < (i)*0.1))
    }
    
    plot(frac,disofrac,type="l",ylim=c(0,1),ylab="%residues",xlab="% disorder")
    lines(frac,Meff100,col="blue")
    lines(frac,Meff200,col="red")
    lines(frac,Meff500,col="green")
    lines(frac,Meff1000,col="grey")

    TM.var=NULL
    TM.var[0]=1
    TM.var[1]=2
    TM.var[2]=3

    TM.PFam=NULL
    TM.PFam[0]=length(which(TM == "M"  & Pfam_ID != "" ))/length(which(TM == "M"))
    TM.PFam[1]=length(which(TM == "i"  & Pfam_ID != "" ))/length(which(TM == "i"))
    TM.PFam[2]=length(which(TM == "o"  & Pfam_ID != "" ))/length(which(TM == "o"))

    TM.Meff100=NULL
    TM.Meff100[0]=length(which(TM == "M"  & Pfam_Meff >100 ))/length(which(TM == "M"))
    TM.Meff100[1]=length(which(TM == "i"  & Pfam_Meff >100  ))/length(which(TM == "i"))
    TM.Meff100[2]=length(which(TM == "o"  & Pfam_Meff >100  ))/length(which(TM == "o"))

    plot(TM.var,TM.Pfam,type="l",ylim=c(0,1),ylab="%residues",xlab="TM",name=c("M","i","o"))
    lines(TM.Meff100,Meff100,col="blue")
    lines(TM.Meff200,Meff200,col="red")
    lines(TM.Meff500,Meff500,col="green")
    lines(TM.Meff1000,Meff1000,col="grey")

    
Pfam_Meff[is.na(Pfam_Meff)] <- 1
Pfam_Meff[Pfam_Meff < 1] <- 1
#Pfam.log <- log10(Pfam_Meff)
    
disoMeff0<-Pfam_Meff[Disorder > 0][Disorder <0.1]
disoMeff1<-Pfam_Meff[Disorder > 0.1][Disorder <0.2]
disoMeff2<-Pfam_Meff[Disorder > 0.2][Disorder <0.3]
disoMeff3<-Pfam_Meff[Disorder > 0.3][Disorder <0.4]
disoMeff4<-Pfam_Meff[Disorder > 0.4][Disorder <0.5]
disoMeff5<-Pfam_Meff[Disorder > 0.5][Disorder <0.6]
disoMeff6<-Pfam_Meff[Disorder > 0.6][Disorder <0.7]
disoMeff7<-Pfam_Meff[Disorder > 0.7][Disorder <0.8]
disoMeff8<-Pfam_Meff[Disorder > 0.8][Disorder <0.9]
disoMeff9<-Pfam_Meff[Disorder > 0.9][Disorder <1.0]
boxplot(disoMeff0,disoMeff1,disoMeff2,disoMeff3,disoMeff4,disoMeff5,disoMeff6,disoMeff7,disoMeff8,disoMeff9,names=c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95),log="y")

diso.1=
diso.2=length(which(Disorder < 0.5) & which(Pfam_ID != ""))/length(which(Disorder > 0.5))

araH.0<-araH$V3[araH$V4>=4.0 ][ araH$V4<4.2 ]

    

detach(dat)
