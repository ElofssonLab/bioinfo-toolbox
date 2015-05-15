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
labels=NULL
loop=seq(1,6)



                                        # Ecoli
ecoli$Pfam_Meff[ecoli$Pfam_Meff==-1 & ecoli$Pfam_pos != -1)] <- 0
ecoli$Pfam_Meff[ecoli$Pfam_Meff==-1 & ecoli$Pfam_pos == -1)] <- 100
ecoli$Pfam_Meff[is.na(ecoli$Pfam_Meff)] <- 0

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
    ecoliPDBNoTM[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$TM == 0  & ecoli$PDB_ID == ""))
    ecoliPDBDiso[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Disorder>=0.5  & ecoli$PDB_ID == ""))
    ecoliPDBNoDiso[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Disorder<0.5  & ecoli$PDB_ID == ""))
    ecoliPDBSeg[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Seg_low > 0 & ecoli$PDB_ID == ""))
    ecoliPDBNoSeg[i]=length(which(ecoli$Pfam_Meff > cutoffs[i] & ecoli$Pfam_Meff <= cutoffs[i+1] & ecoli$Seg_low == 0  & ecoli$PDB_ID == ""))
}

ecoliPDBAll[7]=length(which(ecoli$PDB_ID != "") )
ecoliDisoAll=c(ecoliPDBNoDiso,ecoliPDBDiso,ecoliPDBAll[7])
ecoliPDBTM[7]=length(which(ecoli$PDB_ID != "" & ecoli$TM > 0) )
ecoliPDBNoTM[7]=length(which(ecoli$PDB_ID != "" & ecoli$TM == 0 ) )
ecoliPDBDiso[7]=length(which(ecoli$PDB_ID != "" & ecoli$Disorder>=0.5 ) )
ecoliPDBNoDiso[7]=length(which(ecoli$PDB_ID != ""  & ecoli$Disorder<0.5 ) )
ecoliPDBSeg[7]=length(which(ecoli$PDB_ID != "" & ecoli$Seg_low > 0) )
ecoliPDBNoSeg[7]=length(which(ecoli$PDB_ID != ""  & ecoli$Seg_low == 0) )

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
outfile=paste(genome,"-PDB-pie.png",sep="")
png(outfile)
fraction=(sum(ecoliPDB)/(sum(ecoliNoPDB)+sum(ecoliPDB)))**2
pie(ecoliAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
floating.pie(0,0,ecoliPDB, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliNoPDB, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"PDB",cex=0.4,border=NA)
legend(-.4,.0,"NoPDB",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-TM-pie.png",sep="")
png(outfile)
pie(ecoliAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
floating.pie(0,0,ecoliTM, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"TM",cex=0.4,border=NA)
legend(-.4,.0,"NoTM",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-Diso-pie.png",sep="")
png(outfile)
pie(ecoliAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
floating.pie(0,0,ecoliDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"Diso",cex=0.4,border=NA)
legend(-.4,.0,"NoDiso",cex=0.4,border=NA)
dev.off()

labels[7]="PDB"
colors[7]="green"

genome="ecoli-all"
pct <- round(ecoliPDBAll/sum(ecoliAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste(genome,"-TM-pie.png",sep="")
png(outfile)
pie(ecoliPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
floating.pie(0,0,ecoliPDBTM, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliPDBNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"TM",cex=0.4,border=NA)
legend(-.4,.0,"NoTM",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-Diso-pie.png",sep="")
png(outfile)
pie(ecoliPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
floating.pie(0,0,ecoliPDBDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliPDBNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"Diso",cex=0.4,border=NA)
legend(-.4,.0,"NoDiso",cex=0.4,border=NA)
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
labels=NULL
loop=seq(1,6)



# Sacch
sacch$Pfam_Meff[sacch$Pfam_Meff==-1 & sacch$Pfam_pos != -1)] <- 0
sacch$Pfam_Meff[sacch$Pfam_Meff==-1 & sacch$Pfam_pos == -1)] <- 100
sacch$Pfam_Meff[is.na(sacch$Pfam_Meff)] <- 0

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
    sacchPDBDiso[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Disorder>=0.5  & sacch$PDB_ID == ""))
    sacchPDBNoDiso[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Disorder<0.5  & sacch$PDB_ID == ""))
    sacchPDBSeg[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Seg_low > 0 & sacch$PDB_ID == ""))
    sacchPDBNoSeg[i]=length(which(sacch$Pfam_Meff > cutoffs[i] & sacch$Pfam_Meff <= cutoffs[i+1] & sacch$Seg_low == 0  & sacch$PDB_ID == ""))
}
sacchPDBAll[7]=length(which(sacch$PDB_ID != "") )
sacchDisoAll=c(sacchPDBNoDiso,sacchPDBDiso,sacchPDBAll[7])
sacchPDBTM[7]=length(which(sacch$PDB_ID != "" & sacch$TM > 0) )
sacchPDBNoTM[7]=length(which(sacch$PDB_ID != "" & sacch$TM == 0 ) )
sacchPDBDiso[7]=length(which(sacch$PDB_ID != "" & sacch$Disorder>=0.5 ) )
sacchPDBNoDiso[7]=length(which(sacch$PDB_ID != ""  & sacch$Disorder<0.5 ) )
sacchPDBSeg[7]=length(which(sacch$PDB_ID != "" & sacch$Seg_low > 0) )
sacchPDBNoSeg[7]=length(which(sacch$PDB_ID != ""  & sacch$Seg_low == 0) )

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
outfile=paste(genome,"-PDB-pie.png",sep="")
png(outfile)
fraction=(sum(sacchPDB)/(sum(sacchNoPDB)+sum(sacchPDB)))**2
pie(sacchAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
floating.pie(0,0,sacchPDB, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchNoPDB, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"PDB",cex=0.4,border=NA)
legend(-.4,.0,"NoPDB",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-TM-pie.png",sep="")
png(outfile)
pie(sacchAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
floating.pie(0,0,sacchTM, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"TM",cex=0.4,border=NA)
legend(-.4,.0,"NoTM",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-Diso-pie.png",sep="")
png(outfile)
pie(sacchAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
floating.pie(0,0,sacchDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"Diso",cex=0.4,border=NA)
legend(-.4,.0,"NoDiso",cex=0.4,border=NA)
dev.off()

labels[7]="PDB"
colors[7]="green"

genome="sacch-all"
pct <- round(sacchAll/sum(sacchPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste(genome,"-TM-pie.png",sep="")
png(outfile)
pie(sacchPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
floating.pie(0,0,sacchPDBTM, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchPDBNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"TM",cex=0.4,border=NA)
legend(-.4,.0,"NoTM",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-Diso-pie.png",sep="")
png(outfile)
pie(sacchPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
floating.pie(0,0,sacchPDBDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,sacchPDBNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"Diso",cex=0.4,border=NA)
legend(-.4,.0,"NoDiso",cex=0.4,border=NA)
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
labels=NULL
loop=seq(1,6)



# Homo
homo$Pfam_Meff[homo$Pfam_Meff==-1 & homo$Pfam_pos != -1)] <- 0
homo$Pfam_Meff[homo$Pfam_Meff==-1 & homo$Pfam_pos == -1)] <- 100
homo$Pfam_Meff[is.na(homo$Pfam_Meff)] <- 0


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
    homoPDBDiso[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Disorder>=0.5  & homo$PDB_ID == ""))
    homoPDBNoDiso[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Disorder<0.5  & homo$PDB_ID == ""))
    homoPDBSeg[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Seg_low > 0 & homo$PDB_ID == ""))
    homoPDBNoSeg[i]=length(which(homo$Pfam_Meff > cutoffs[i] & homo$Pfam_Meff <= cutoffs[i+1] & homo$Seg_low == 0  & homo$PDB_ID == ""))
}
homoPDBAll[7]=length(which(homo$PDB_ID != "") )
homoDisoAll=c(homoPDBNoDiso,homoPDBDiso,homoPDBAll[7])
homoPDBTM[7]=length(which(homo$PDB_ID != "" & homo$TM > 0) )
homoPDBNoTM[7]=length(which(homo$PDB_ID != "" & homo$TM == 0 ) )
homoPDBDiso[7]=length(which(homo$PDB_ID != "" & homo$Disorder>=0.5 ) )
homoPDBNoDiso[7]=length(which(homo$PDB_ID != ""  & homo$Disorder<0.5 ) )
homoPDBSeg[7]=length(which(homo$PDB_ID != "" & homo$Seg_low > 0) )
homoPDBNoSeg[7]=length(which(homo$PDB_ID != ""  & homo$Seg_low == 0) )

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
outfile=paste(genome,"-PDB-pie.png",sep="")
png(outfile)
fraction=(sum(homoPDB)/(sum(homoNoPDB)+sum(homoPDB)))**2
pie(homoAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
floating.pie(0,0,homoPDB, col=colors,main='',radius=0.7)
floating.pie(0,0,homoNoPDB, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"PDB",cex=0.4,border=NA)
legend(-.4,.0,"NoPDB",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-TM-pie.png",sep="")
png(outfile)
pie(homoAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
floating.pie(0,0,homoTM, col=colors,main='',radius=0.7)
floating.pie(0,0,homoNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"TM",cex=0.4,border=NA)
legend(-.4,.0,"NoTM",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-Diso-pie.png",sep="")
png(outfile)
pie(homoAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
floating.pie(0,0,homoDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,homoNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"Diso",cex=0.4,border=NA)
legend(-.4,.0,"NoDiso",cex=0.4,border=NA)
dev.off()

labels[7]="PDB"
colors[7]="green"

genome="homo-all"
pct <- round(homoAll/sum(homoPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste(genome,"-TM-pie.png",sep="")
png(outfile)
pie(homoPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
floating.pie(0,0,homoPDBTM, col=colors,main='',radius=0.7)
floating.pie(0,0,homoPDBNoTM, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"TM",cex=0.4,border=NA)
legend(-.4,.0,"NoTM",cex=0.4,border=NA)
dev.off()

outfile=paste(genome,"-Diso-pie.png",sep="")
png(outfile)
pie(homoPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
floating.pie(0,0,homoPDBDiso, col=colors,main='',radius=0.7)
floating.pie(0,0,homoPDBNoDiso, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
#legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
legend(-1,.0,"all",cex=0.4,border=NA)
legend(-.7,.0,"Diso",cex=0.4,border=NA)
legend(-.4,.0,"NoDiso",cex=0.4,border=NA)
dev.off()



                                        #Summary of all
genome="All"
outfile=paste(genome,"-pie.png",sep="")
png(outfile)

pie(homoPDBAll, labels=labels,col=colors,main=genome,radius=iniR,border = NA)
floating.pie(0,0,sacchPDBAll, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliPDBAll, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"Homo",cex=0.4,border=NA)
legend(-.7,.0,"Sacc",cex=0.4,border=NA)
legend(-.4,.0,"EColi",cex=0.4,border=NA)
dev.off()


labels[7]="Diso-NoPfam"
colors[7]="lightpink"
labels[8]="Pfam diso"
colors[8]="lightpink2"
labels[9]=">100 diso"
colors[9]="deeppink1"
labels[10]=">200 diso"
colors[10]="deeppink2"
labels[11]=">500 diso"
colors[11]="deeppink3"
labels[12]=">1000 diso"
colors[12]="red"
labels[13]="PDB"
colors[13]="green"

#ecoliDisoAll=c(ecoliPDBNoDiso,ecoliPDBDiso,ecoliPDBAll[7])
#sacchDisoAll=c(sacchPDBNoDiso,sacchPDBDiso,sacchPDBAll[7])
#homoDisoAll=c(homoPDBNoDiso,homoPDBDiso,homoPDBAll[7])

genome="All"
outfile=paste(genome,"Diso-pie.png",sep="")

pct <- round(ecoliDisoAll/sum(ecoliDisoAll)*100,digits=1)
Names <- paste(labels,pct)
names <- paste(names,"%",sep="")
png(outfile)

pie(homoDisoAll, labels=labels,col=colors,main=genome,radius=iniR,border = NA)
floating.pie(0,0,sacchDisoAll, col=colors,main='',radius=0.7)
floating.pie(0,0,ecoliDisoAll, col=colors,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"Homo",cex=0.4,border=NA)
legend(-.7,.0,"Sacch",cex=0.4,border=NA)
legend(-.4,.0,"EColi",cex=0.4,border=NA)
dev.off()
