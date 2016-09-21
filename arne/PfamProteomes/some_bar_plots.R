library(vioplot)
library('plotrix')


df=NULL
dfcov=NULL
test=NULL
genomes=NULL
#genomes<-append(genomes,"escherichia_coli.df")
#genomes<-append(genomes,"saccharomyces_cerevisae.df")
#genomes<-append(genomes,"homo_sapiens.df")
#genomes<-append(genomes,"escherichia_coli.df.2015.blast.29.0")
genomes<-append(genomes,"escherichia_coli_ehec.df.2015.blast.29.0")
genomes<-append(genomes,"helicobacter_pylori.df.2015.blast.29.0")
genomes<-append(genomes,"staphylococcus_aureus_mrsa.df.2015.blast.29.0")
num=0
                                        #genome="escherichia_coli"
for (genome in genomes){
num=num+1
file=paste("data/",genome,".tsv",sep="")

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


                                        # We need to do this in a better way !!

                                        # First make some new columns


#------------------------------------------------------------------------

                                        #plot(Disorder,dat$Pfam_Meff)
#Pfam=dat$Pfam_Meff(replaace

#Fraction disorder with PFam hit
loop=seq(1,10)
disoPfam=NULL
disofrac=NULL
Meff100=NULL
Meff1000=NULL
frac=NULL
for (i in loop){
    frac[i]=0.1*i
    disofrac[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))/length(dat$Disorder)
    disoPfam[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1  & dat$Pfam_ID != ""))/length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))
    Meff100[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1  & dat$Pfam_Meff >100 ))/length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))
    Meff1000[i]=length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1  & dat$Pfam_Meff >1000 ))/length(which(dat$Disorder > (i-1)*0.1 & dat$Disorder < (i)*0.1))
}


ntrees <- 10
linetype <- c(1:ntrees)




#outfile=paste(genome,"-DisoBox.pdf",sep="")
#pdf(outfile)
#boxplot(disoMeff0,disoMeff1,disoMeff2,disoMeff3,disoMeff4,disoMeff5,disoMeff6,disoMeff7,disoMeff8,disoMeff9,names=c("0.05","0.15","0.25","0.35","0.45","0.55","0.65","0.75","0.85","0.95"),log="y",main=genome,ylab="%residues",xlab="% disorder")
#dev.off()


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

TMMeff1000=NULL
TMMeff1000[1]=length(which(dat$TM_topogy== "M"  & dat$Pfam_Meff >1000 ))/length(which(dat$TM_topogy== "M"))
TMMeff1000[2]=length(which(dat$TM_topogy== "i"  & dat$Pfam_Meff >1000 ))/length(which(dat$TM_topogy== "i"))
TMMeff1000[3]=length(which(dat$TM_topogy== "o"  & dat$Pfam_Meff >1000 ))/length(which(dat$TM_topogy== "o"))



# Pie Charts

# Fraction of residues with Pfam-coverage and Meff coverage=NULL
coverage=NULL
labels=NULL
coverage[1]=length(which(dat$Pfam_ID == "" ))
labels[1]="NoPfam"
coverage[2]=length(which(dat$Pfam_Meff > 1000))
labels[2]=">1000"
coverage[3]=length(which(dat$Pfam_Meff > 100 & dat$Pfam_Meff <= 1000))
labels[3]=">100"
coverage[4]=length(which(dat$Pfam_ID != "" & dat$Pfam_Meff<=100))
labels[4]="<100"
pct <- round(coverage/sum(coverage)*100)
outfile=paste(genome,"-coverage-new.pdf",sep="")
pdf(outfile)
pie(coverage,labels = labels, col=rainbow(length(labels)),main=genome)
dev.off()





# Now for PDB info 

# Pie Charts

# Fraction of residues with Pfam-coverage and Meff coverage=NULL
coverage=NULL
labels=NULL
colors=NULL
coverage<-append(coverage,length(which(dat$PDB_ID != "" & dat$PDB_E.value == 0)))
labels<-append(labels,"PDB")
colors<-append(colors,"black")

coverage<-append(coverage,length(which(dat$PDB_ID != "" & dat$PDB_E.value > 0 )))
labels<-append(labels,"Homology")
colors<-append(colors,"grey80")

coverage<-append(coverage,length(which(dat$Pfam_Meff > 1000 & dat$PDB_ID == "")))
labels<-append(labels,"Meff >1000")
colors<-append(colors,"Red")

coverage<-append(coverage,length(which(dat$Pfam_Meff > 100 & dat$Pfam_Meff <= 1000 & dat$PDB_ID == "")))
labels<-append(labels,"Meff >100")
colors<-append(colors,"green")

coverage<-append(coverage,length(which(dat$Pfam_ID != "" & dat$Pfam_Meff <= 100 & dat$PDB_ID == "")))
labels<-append(labels,"Meff < 100")
colors<-append(colors,"lightblue")


coverage<-append(coverage,length(which(dat$Pfam_ID == "" & dat$PDB_ID == "" )))
labels<-append(labels,"NoPfam")
colors<-append(colors,"white")

pct <- (coverage/sum(coverage))
outfile=paste(genome,"-coverage-PDB2.png",sep="")
png(outfile)
pie(coverage,labels = labels, col=rainbow(length(labels)),main=genome)
dev.off()

df = append(df,pct)
dfcov = append(dfcov,coverage)


}    

                                        #Pfam.log <- log10(dat$Pfam_Meff)
rows=length(colors)
test=matrix(c(df*100),nrow=rows,ncol=num)

outfile=paste("figures/genomes-bar.png",sep="")
png(outfile,width=1280,height=1280)
barplot(test,main="Precentage of Residues",col=colors,legend=(labels),xlim=c(0,c(num+.5)),xlab="Specie",ylab="Percentage of residues",cex.names=1,cex.axis=1.,cex=1,names=genomes)
dev.off()
