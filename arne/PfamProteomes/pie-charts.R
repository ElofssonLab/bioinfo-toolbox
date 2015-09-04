library(vioplot)
library('plotrix')
genomes=NULL
genomes[1]="escherichia_coli"
genomes[2]="saccharomyces_cerevisae"
genomes[3]="homo_sapiens"
#genome="escherichia_coli"
#for (genome in genomes){
file=paste("data/",genomes[1],".df.tsv",sep="")
ecoli<-read.table(file, sep='\t', header=T)
file=paste("data/",genomes[2],".df.tsv",sep="")
sacch<-read.table(file, sep='\t', header=T)
file=paste("data/",genomes[3],".df.tsv",sep="")
homo<-read.table(file, sep='\t', header=T)
#}

columns=NULL
cutoffs=NULL

cutoffs[1]=-2
cutoffs[2]=0
cutoffs[3]=99
cutoffs[4]=199
cutoffs[5]=499
cutoffs[6]=999
cutoffs[7]=100000000000

# --- Preprocessing --
                                        # Ecoli
ecoli$Pfam_Meff[ecoli$Pfam_Meff==-1 & ecoli$Pfam_pos != -1] <- 0
ecoli$Pfam_Meff[ecoli$Pfam_Meff==-1 & ecoli$Pfam_pos == -1] <- 100
ecoli$Pfam_Meff[is.na(ecoli$Pfam_Meff)] <- 0

# Sacch
sacch$Pfam_Meff[sacch$Pfam_Meff==-1 & sacch$Pfam_pos != -1] <- 0
sacch$Pfam_Meff[sacch$Pfam_Meff==-1 & sacch$Pfam_pos == -1] <- 100
sacch$Pfam_Meff[is.na(sacch$Pfam_Meff)] <- 0

# Homo
homo$Pfam_Meff[homo$Pfam_Meff==-1 & homo$Pfam_pos != -1] <- 0
homo$Pfam_Meff[homo$Pfam_Meff==-1 & homo$Pfam_pos == -1] <- 100
homo$Pfam_Meff[is.na(homo$Pfam_Meff)] <- 0


#-----------------------  E. Coli  --------------------------------
ecoliAll=NULL
ecoliTM=NULL
ecoliDiso=NULL
ecoliPDB=NULL
ecoliSeg=NULL
ecoliNoTM=NULL
ecoliNoDiso=NULL
ecoliNoPDB=NULL
ecoliNoSeg=NULL
ecoliPDBAll=NULL
ecoliPDBTM=NULL
ecoliPDBDiso=NULL
ecoliPDBPDB=NULL
ecoliPDBSeg=NULL
ecoliPDBNoTM=NULL
ecoliPDBNoDiso=NULL
ecoliPDBNoSeg=NULL
ecoliPDBTMNoDiso=NULL
ecoliPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)




                                        # Ja
for (i in loop){
    ecoliAll[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1]))
    ecoliTM[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$TM > 0))
    ecoliNoTM[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$TM == 0))
    ecoliDiso[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Disorder>=0.5))
    ecoliNoDiso[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Disorder<0.5))
    ecoliPDB[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$PDB_ID != ""))
    ecoliNoPDB[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$PDB_ID == ""))
    ecoliSeg[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Seg_low > 0))
    ecoliNoSeg[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Seg_low == 0))
    labels[i]=cutoffs[i]
    ecoliPDBAll[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$PDB_ID == "") )
    ecoliPDBTM[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$TM > 0 & ecoli$PDB_ID == ""))
    ecoliPDBNoTM[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$TM == 0  & ecoli$PDB_ID == "" & ecoli$Disorder<0.5 ))
    ecoliPDBTMNoDiso[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$TM > 0 & ecoli$PDB_ID == ""& ecoli$Disorder<0.5 ))
    ecoliPDBNoTMNoDiso[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$TM == 0  & ecoli$PDB_ID == ""))
    ecoliPDBDiso[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Disorder>=0.5  & ecoli$PDB_ID == ""))
    ecoliPDBNoDiso[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Disorder<0.5  & ecoli$PDB_ID == ""))
    ecoliPDBSeg[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Seg_low > 0 & ecoli$PDB_ID == ""))
    ecoliPDBNoSeg[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Seg_low == 0  & ecoli$PDB_ID == ""))
}

ecoliPDBAll[7]=length(which(ecoli$PDB_ID != "") )
ecoliDisoAll=c(ecoliPDBNoDiso,ecoliPDBDiso,ecoliPDBAll[7])
ecoliTMAll=c(ecoliPDBNoTMNoDiso[1:6],ecoliPDBTMNoDiso[1:6],ecoliPDBDiso[1:6],ecoliPDB[1:6])
ecoliPDBTM[7]=length(which(ecoli$PDB_ID != "" & ecoli$TM > 0) )
ecoliPDBNoTM[7]=length(which(ecoli$PDB_ID != "" & ecoli$TM == 0 ) )
ecoliPDBDiso[7]=length(which(ecoli$PDB_ID != "" & ecoli$Disorder>=0.5 ) )
ecoliPDBNoDiso[7]=length(which(ecoli$PDB_ID != ""  & ecoli$Disorder<0.5 ) )
ecoliPDBSeg[7]=length(which(ecoli$PDB_ID != "" & ecoli$Seg_low > 0) )
ecoliPDBNoSeg[7]=length(which(ecoli$PDB_ID != ""  & ecoli$Seg_low == 0) )


ecoliAll[ecoliAll==0]<-0.001
ecoliPDBAll[ecoliPDBAll==0]<-0.001
ecoliDisoAll[ecoliDisoAll==0]<-0.001
ecoliTMAll[ecoliTMAll==0]<-0.001

labels=NULL
colors=NULL
labels[1]="NoPfam"
colors[1]="cadetblue1"
labels[2]="Pfam"
colors[2]="cadetblue2"
labels[3]=">100"
colors[3]="blue1"
labels[4]=">200"
colors[4]="blue2"
labels[5]=">500"
colors[5]="blue3"
labels[6]=">1000"
colors[6]="darkblue"


pct <- round(ecoliAll/sum(ecoliAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="ecoli"
outfile=paste("figures/",genome,"-PDB-pie.png",sep="")
png(outfile,width=1280,height=1280)
fraction=(sum(ecoliPDB)/(sum(ecoliNoPDB)+sum(ecoliPDB)))**2
pie(ecoliAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,ecoliPDB, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliNoPDB, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"PDB",cex=2.0,border=NA)
legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-TM-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(ecoliAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,ecoliTM, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"TM",cex=2.0,border=NA)
legend(-.4,.0,"NoTM",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Diso-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(ecoliAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,ecoliDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"Diso",cex=2.0,border=NA)
legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
dev.off()

labels[7]="PDB"
colors[7]="green"

genome="ecoli-all"
pct <- round(ecoliPDBAll/sum(ecoliAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,"-TM-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(ecoliPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,ecoliPDBTM, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliPDBNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"TM",cex=2.0,border=NA)
legend(-.4,.0,"NoTM",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Diso-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(ecoliPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,ecoliPDBDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliPDBNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"Diso",cex=2.0,border=NA)
legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
dev.off()


#-----------------------  Sacharomyces  --------------------------------
sacchAll=NULL
sacchTM=NULL
sacchDiso=NULL
sacchPDB=NULL
sacchSeg=NULL
sacchNoTM=NULL
sacchNoDiso=NULL
sacchNoPDB=NULL
sacchNoSeg=NULL
sacchPDBAll=NULL
sacchPDBTM=NULL
sacchPDBDiso=NULL
sacchPDBPDB=NULL
sacchPDBSeg=NULL
sacchPDBNoTM=NULL
sacchPDBNoDiso=NULL
sacchPDBNoSeg=NULL
sacchPDBTMNoDiso=NULL
sacchPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)




                                        # Ja
for (i in loop){
    sacchAll[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1]))
    sacchTM[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$TM > 0))
    sacchNoTM[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$TM == 0))
    sacchDiso[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Disorder>=0.5))
    sacchNoDiso[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Disorder<0.5))
    sacchPDB[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$PDB_ID != ""))
    sacchNoPDB[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$PDB_ID == ""))
    sacchSeg[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Seg_low > 0))
    sacchNoSeg[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Seg_low == 0))
    labels[i]=cutoffs[i]
    sacchPDBAll[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$PDB_ID == "") )
    sacchPDBTM[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$TM > 0 & sacch$PDB_ID == ""))
    sacchPDBNoTM[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$TM == 0  & sacch$PDB_ID == ""))
    sacchPDBTMNoDiso[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$TM > 0 & sacch$PDB_ID == ""& sacch$Disorder<0.5 ))
    sacchPDBNoTMNoDiso[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$TM == 0  & sacch$PDB_ID == "" & sacch$Disorder<0.5 ))
    sacchPDBDiso[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Disorder>=0.5  & sacch$PDB_ID == ""))
    sacchPDBNoDiso[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Disorder<0.5  & sacch$PDB_ID == ""))
    sacchPDBSeg[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Seg_low > 0 & sacch$PDB_ID == ""))
    sacchPDBNoSeg[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Seg_low == 0  & sacch$PDB_ID == ""))
}
sacchPDBAll[7]=length(which(sacch$PDB_ID != "") )
sacchTMAll=c(sacchPDBNoTMNoDiso[1:6],sacchPDBTMNoDiso[1:6],sacchPDBDiso[1:6],sacchPDB[1:6])
sacchDisoAll=c(sacchPDBNoDiso,sacchPDBDiso,sacchPDBAll[7])
sacchPDBTM[7]=length(which(sacch$PDB_ID != "" & sacch$TM > 0) )
sacchPDBNoTM[7]=length(which(sacch$PDB_ID != "" & sacch$TM == 0 ) )
sacchPDBDiso[7]=length(which(sacch$PDB_ID != "" & sacch$Disorder>=0.5 ) )
sacchPDBNoDiso[7]=length(which(sacch$PDB_ID != ""  & sacch$Disorder<0.5 ) )
sacchPDBSeg[7]=length(which(sacch$PDB_ID != "" & sacch$Seg_low > 0) )
sacchPDBNoSeg[7]=length(which(sacch$PDB_ID != ""  & sacch$Seg_low == 0) )


sacchAll[sacchAll==0]<-0.001
sacchPDBAll[sacchPDBAll==0]<-0.001
sacchDisoAll[sacchDisoAll==0]<-0.001
sacchTMAll[sacchTMAll==0]<-0.001

labels=NULL
colors=NULL
labels[1]="NoPfam"
colors[1]="cadetblue1"
labels[2]="Pfam"
colors[2]="cadetblue2"
labels[3]=">100"
colors[3]="blue1"
labels[4]=">200"
colors[4]="blue2"
labels[5]=">500"
colors[5]="blue3"
labels[6]=">1000"
colors[6]="darkblue"


pct <- round(sacchAll/sum(sacchAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="sacch"
outfile=paste("figures/",genome,"-PDB-pie.png",sep="")
png(outfile,width=1280,height=1280)
fraction=(sum(sacchPDB)/(sum(sacchNoPDB)+sum(sacchPDB)))**2
pie(sacchAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchPDB, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchNoPDB, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"PDB",cex=2.0,border=NA)
legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-TM-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(sacchAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchTM, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"TM",cex=2.0,border=NA)
legend(-.4,.0,"NoTM",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Diso-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(sacchAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"Diso",cex=2.0,border=NA)
legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
dev.off()

labels[7]="PDB"
colors[7]="green"

genome="sacch-all"
pct <- round(sacchAll/sum(sacchPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,"-TM-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(sacchPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchPDBTM, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchPDBNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"TM",cex=2.0,border=NA)
legend(-.4,.0,"NoTM",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Diso-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(sacchPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchPDBDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchPDBNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"Diso",cex=2.0,border=NA)
legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
dev.off()



#-----------------------  Homo Sapiens  --------------------------------
homoAll=NULL
homoTM=NULL
homoDiso=NULL
homoPDB=NULL
homoSeg=NULL
homoNoTM=NULL
homoNoDiso=NULL
homoNoPDB=NULL
homoNoSeg=NULL
homoPDBAll=NULL
homoPDBTM=NULL
homoPDBDiso=NULL
homoPDBPDB=NULL
homoPDBSeg=NULL
homoPDBNoTM=NULL
homoPDBNoDiso=NULL
homoPDBNoSeg=NULL
homoPDBTMNoDiso=NULL
homoPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)





                                        # Ja
for (i in loop){
    homoAll[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1]))
    homoTM[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$TM > 0))
    homoNoTM[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$TM == 0))
    homoDiso[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Disorder>=0.5))
    homoNoDiso[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Disorder<0.5))
    homoPDB[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$PDB_ID != ""))
    homoNoPDB[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$PDB_ID == ""))
    homoSeg[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Seg_low > 0))
    homoNoSeg[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Seg_low == 0))
    labels[i]=cutoffs[i]
    homoPDBAll[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$PDB_ID == "") )
    homoPDBTM[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$TM > 0 & homo$PDB_ID == ""))
    homoPDBNoTM[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$TM == 0  & homo$PDB_ID == ""))
    homoPDBNoTMNoDiso[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$TM == 0  & homo$PDB_ID == "" & homo$Disorder<0.5 ))
    homoPDBTMNoDiso[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$TM > 0 & homo$PDB_ID == ""& homo$Disorder<0.5 ))
    homoPDBDiso[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Disorder>=0.5  & homo$PDB_ID == ""))
    homoPDBNoDiso[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Disorder<0.5  & homo$PDB_ID == ""))
    homoPDBSeg[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Seg_low > 0 & homo$PDB_ID == ""))
    homoPDBNoSeg[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Seg_low == 0  & homo$PDB_ID == ""))
}
homoPDBAll[7]=length(which(homo$PDB_ID != "") )
homoDisoAll=c(homoPDBNoDiso,homoPDBDiso,homoPDBAll[7])
homoTMAll=c(homoPDBNoTMNoDiso[1:6],homoPDBTMNoDiso[1:6],homoPDBDiso[1:6],homoPDB[1:6])
homoPDBTM[7]=length(which(homo$PDB_ID != "" & homo$TM > 0) )
homoPDBNoTM[7]=length(which(homo$PDB_ID != "" & homo$TM == 0 ) )
homoPDBDiso[7]=length(which(homo$PDB_ID != "" & homo$Disorder>=0.5 ) )
homoPDBNoDiso[7]=length(which(homo$PDB_ID != ""  & homo$Disorder<0.5 ) )
homoPDBSeg[7]=length(which(homo$PDB_ID != "" & homo$Seg_low > 0) )
homoPDBNoSeg[7]=length(which(homo$PDB_ID != ""  & homo$Seg_low == 0) )

homoAll[homoAll==0]<-0.001
homoPDBAll[homoPDBAll==0]<-0.001
homoDisoAll[homoDisoAll==0]<-0.001
homoTMAll[homoTMAll==0]<-0.001

labels=NULL
colors=NULL
labels[1]="NoPfam"
colors[1]="cadetblue1"
labels[2]="Pfam"
colors[2]="cadetblue2"
labels[3]=">100"
colors[3]="blue1"
labels[4]=">200"
colors[4]="blue2"
labels[5]=">500"
colors[5]="blue3"
labels[6]=">1000"
colors[6]="darkblue"


pct <- round(homoAll/sum(homoAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="homo"
outfile=paste("figures/",genome,"-PDB-pie.png",sep="")
png(outfile,width=1280,height=1280)
fraction=(sum(homoPDB)/(sum(homoNoPDB)+sum(homoPDB)))**2
pie(homoAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,homoPDB, col=colors,main='',radius=0.7)
floating.pie(0,0,homoNoPDB, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"PDB",cex=2.0,border=NA)
legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-TM-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(homoAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,homoTM, col=colors,main='',radius=0.7)
floating.pie(0,0,homoNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"TM",cex=2.0,border=NA)
legend(-.4,.0,"NoTM",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Diso-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(homoAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,homoDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,homoNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"Diso",cex=2.0,border=NA)
legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
dev.off()

labels[7]="PDB"
colors[7]="green"

genome="homo-all"
pct <- round(homoAll/sum(homoPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,"-TM-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(homoPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,homoPDBTM, col=colors,main='',radius=0.7)
floating.pie(0,0,homoPDBNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"TM",cex=2.0,border=NA)
legend(-.4,.0,"NoTM",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Diso-pie.png",sep="")
png(outfile,width=1280,height=1280)
pie(homoPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,homoPDBDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,homoPDBNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=2.0,border=NA)
legend(-.7,.0,"Diso",cex=2.0,border=NA)
legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
dev.off()



                                        #Summary of all
genome="All"
outfile=paste("figures/",genome,"-pie.png",sep="")
png(outfile,width=1280,height=1280)

pie(homoPDBAll, labels=labels,col=colors,main=genome,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchPDBAll, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliPDBAll, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"Homo",cex=2.0,border=NA)
legend(-.7,.0,"Sacc",cex=2.0,border=NA)
legend(-.4,.0,"EColi",cex=2.0,border=NA)
dev.off()


labels[7]="TM-NoPfam"
colors[7]="lightpink"
labels[8]="Pfam TM"
colors[8]="lightpink2"
labels[9]=">100 TM"
colors[9]="deeppink1"
labels[10]=">200 TM"
colors[10]="deeppink2"
labels[11]=">500 TM"
colors[11]="deeppink3"
labels[12]=">1000 TM"
colors[12]="deeppink4"
labels[13]="TM-NoPfam"
colors[13]="greenyellow"
labels[14]="Pfam diso"
colors[14]="green"
labels[15]=">100 diso"
colors[15]="green1"
labels[16]=">200 diso"
colors[16]="green3"
labels[17]=">500 diso"
colors[17]="green3"
labels[18]=">1000 diso"
colors[18]="green4"

#ecoliDisoAll=c(ecoliPDBNoDiso,ecoliPDBDiso,ecoliPDBAll[7])
#sacchDisoAll=c(sacchPDBNoDiso,sacchPDBDiso,sacchPDBAll[7])
#homoDisoAll=c(homoPDBNoDiso,homoPDBDiso,homoPDBAll[7])

genome="All"
outfile=paste("figures/",genome,"Diso-pie.png",sep="")

pct <- round(ecoliDisoAll/sum(ecoliDisoAll)*100,digits=1)
Names <- paste(labels,pct)
names <- paste(Names,"%",sep="")
png(outfile,width=1280,height=1280)

pie(homoDisoAll, labels=labels,col=colors,main=genome,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchDisoAll, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliDisoAll, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"Homo",cex=2.0,border=NA)
legend(-.7,.0,"Sacch",cex=2.0,border=NA)
legend(-.4,.0,"EColi",cex=2.0,border=NA)
dev.off()

# WHat is the interesting part.


labels[7]="TM-NoPfam"
colors[7]="lightpink"
labels[8]="Pfam TM"
colors[8]="lightpink2"
labels[9]=">100 TM"
colors[9]="deeppink1"
labels[10]=">200 TM"
colors[10]="deeppink2"
labels[11]=">500 TM"
colors[11]="deeppink3"
labels[12]=">1000 TM"
colors[12]="red"
labels[13]="TM-NoPfam"
colors[13]="greenyellow"
labels[14]="Pfam diso"
colors[14]="green"
labels[15]=">100 diso"
colors[15]="green1"
labels[16]=">200 diso"
colors[16]="green3"
labels[17]=">500 diso"
colors[17]="green3"
labels[18]=">1000 diso"
colors[18]="green4"
labels[19]="PDB NoPfam"
colors[19]="grey90"
labels[20]="PDB Pfam"
colors[20]="grey80"
labels[21]="PDB >100 "
colors[21]="grey40"
colors[22]="grey30"
labels[22]="PDB >200"
colors[23]="grey20"
labels[23]="PDB >500"
colors[24]="black"
labels[24]="PDB >1000"


Labels=NULL

Labels[1]="Globular"
Labels[7]="TM"
Labels[13]="Disorder"
Labels[19]="PDB"


genome="All"
outfile=paste("figures/",genome,"Diso-pie.png",sep="")
pct <- round(ecoliTMAll/sum(ecoliTMAll)*100,digits=1)
Names <- paste(labels,pct)
names <- paste(Names,"%",sep="")
png(outfile,width=1280,height=1280)
pie(homoTMAll, labels=labels,col=colors,main=genome,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchTMAll, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliTMAll, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"Homo",cex=2.0,border=NA)
legend(-.7,.0,"Sacch",cex=2.0,border=NA)
legend(-.4,.0,"EColi",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-bar.png",sep="")
png(outfile,width=1280,height=1280)
homofrac<-homoTMAll/sum(homoTMAll)
sacchfrac<-sacchTMAll/sum(sacchTMAll)
ecolifrac<-ecoliTMAll/sum(ecoliTMAll)
genomes=NULL
genomes[1]="Homo sapiens"
genomes[2]="Yeast"
genomes[3]="E. Coli"
test=matrix(c(rev(homofrac),rev(sacchfrac),rev(ecolifrac)),nrow=24,ncol=3)
bp <-barplot(test,col=rev(colors),main="Fraction of residues",legend=rev(names),xlim=c(0,4.5),xlab="",ylab="Fraction of residues",names=genomes,cex.names=3,cex.axis=3.,cex=3)


                                        # Calculate fractions
EcoliNoModel=sum(ecoliTMAll[1:2],ecoliTMAll[7:8],ecoliTMAll[13:14])/sum(ecoliTMAll)*100
EcoliPDB=sum(ecoliTMAll[19:24])/sum(ecoliTMAll)*100
EcoliGLOB=sum(ecoliTMAll[3:6])/sum(ecoliTMAll)*100
EcoliTM=sum(ecoliTMAll[9:12])/sum(ecoliTMAll)*100
EcoliTMall=sum(ecoliTMAll[7:12])/sum(ecoliTMAll)*100
EcoliDISO=sum(ecoliTMAll[15:18])/sum(ecoliTMAll)*100
EcoliDISOall=sum(ecoliTMAll[13:18])/sum(ecoliTMAll)*100
EcoliHundred=sum(ecoliTMAll[3:6],ecoliTMAll[9:12],ecoliTMAll[15:18])/sum(ecoliTMAll)*100
EcoliThousand=sum(ecoliTMAll[6],ecoliTMAll[12],ecoliTMAll[18])/sum(ecoliTMAll)*100
EcoliPDB
EcoliHundred
EcoliPDB+EcoliHundred
                                        # Calculate fractions
SacchNoModel=sum(sacchTMAll[1:2],sacchTMAll[7:8],sacchTMAll[13:14])/sum(sacchTMAll)*100
SacchPDB=sum(sacchTMAll[19:24])/sum(sacchTMAll)*100
SacchGLOB=sum(sacchTMAll[3:6])/sum(sacchTMAll)*100
SacchTM=sum(sacchTMAll[9:12])/sum(sacchTMAll)*100
SacchTMall=sum(sacchTMAll[7:12])/sum(sacchTMAll)*100
SacchDISO=sum(sacchTMAll[15:18])/sum(sacchTMAll)*100
SacchDISOall=sum(sacchTMAll[12:18])/sum(sacchTMAll)*100
SacchHundred=sum(sacchTMAll[3:6],sacchTMAll[9:12],sacchTMAll[15:18])/sum(sacchTMAll)*100
SacchThousand=sum(sacchTMAll[6],sacchTMAll[12],sacchTMAll[18])/sum(sacchTMAll)*100
SacchPDB
SacchHundred
SacchPDB+SacchHundred
                                        # Calculate fractions
HomoNoModel=sum(homoTMAll[1:2],homoTMAll[7:8],homoTMAll[13:14])/sum(homoTMAll)*100
HomoPDB=sum(homoTMAll[19:24])/sum(homoTMAll)*100
HomoGLOB=sum(homoTMAll[3:6])/sum(homoTMAll)*100
HomoTM=sum(homoTMAll[9:12])/sum(homoTMAll)*100
HomoTMall=sum(homoTMAll[7:12])/sum(homoTMAll)*100
HomoDISO=sum(homoTMAll[15:18])/sum(homoTMAll)*100
HomoDISOall=sum(homoTMAll[13:18])/sum(homoTMAll)*100
HomoHundred=sum(homoTMAll[3:6],homoTMAll[9:12],homoTMAll[15:18])/sum(homoTMAll)*100
HomoThousand=sum(homoTMAll[6],homoTMAll[12],homoTMAll[18])/sum(homoTMAll)*100
HomoPDB
HomoHundred
HomoPDB+HomoHundred

fractionPDB=NULL
fractionPDB[1]=round(HomoPDB,0)
fractionPDB[2]=round(SacchPDB,0)
fractionPDB[3]=round(EcoliPDB,0)
labelPDB <- paste(fractionPDB,"%",sep="")
text(bp, 0, labelPDB,cex=2,pos=3,col="white") 



fractionDISO=NULL
fractionDISOall=NULL
fractionDISO[1]=round(HomoDISO,1)
fractionDISO[2]=round(SacchDISO,1)
fractionDISO[3]=round(EcoliDISO,1)
fractionDISOall[1]=round(HomoDISOall,1)
fractionDISOall[2]=round(SacchDISOall,1)
fractionDISOall[3]=round(EcoliDISOall,1)
labelDISO <- paste(fractionDISO,"%",sep="")
text(bp, fractionPDB/100 , labelDISO,cex=2,pos=3,col="black") 


fractionTM=NULL
fractionTMall=NULL
fractionTM[1]=round(HomoTM,0)
fractionTM[2]=round(SacchTM,0)
fractionTM[3]=round(EcoliTM,0)
fractionTMall[1]=round(HomoTMall,1)
fractionTMall[2]=round(SacchTMall,1)
fractionTMall[3]=round(EcoliTMall,1)
labelTM <- paste(fractionTM,"%",sep="")
text(bp, (fractionPDB+fractionDISOall )/100, labelTM,cex=2,pos=3,col="black") 

fractionGLOB=NULL
fractionGLOBall=NULL
fractionGLOB[1]=round(HomoGLOB,0)
fractionGLOB[2]=round(SacchGLOB,0)
fractionGLOB[3]=round(EcoliGLOB,0)
labelGLOB <- paste(fractionGLOB,"%",sep="")
text(bp, ( fractionPDB+fractionDISOall+fractionTMall )/100,labelGLOB,cex=2,pos=3,col="red") 

dev.off()
