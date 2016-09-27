library('vioplot')
library('plotrix')
genomes=NULL
years=NULL
years[1]="2005"
years[2]="2010"
years[3]="2015"
genomes[1]="escherichia_coli"
genomes[2]="saccharomyces_cerevisae"
genomes[3]="homo_sapiens"

year=years[1]
file=paste("data/",genomes[1],".df.",year,".blast.prot.tsv",sep="")
ecoliA<-read.table(file, sep='\t', header=T)
file=paste("data/",genomes[2],".df.",year,".blast.prot.tsv",sep="")
sacchA<-read.table(file, sep='\t', header=T)
file=paste("data/",genomes[3],".df.",year,".blast.prot.tsv",sep="")
homoA<-read.table(file, sep='\t', header=T)

year=years[2]
file=paste("data/",genomes[1],".df.",year,".blast.prot.tsv",sep="")
ecoliB<-read.table(file, sep='\t', header=T)
file=paste("data/",genomes[2],".df.",year,".blast.prot.tsv",sep="")
sacchB<-read.table(file, sep='\t', header=T)
file=paste("data/",genomes[3],".df.",year,".blast.prot.tsv",sep="")
homoB<-read.table(file, sep='\t', header=T)

year=years[3]
file=paste("data/",genomes[1],".df.",year,".blast.prot.tsv",sep="")
ecoliC<-read.table(file, sep='\t', header=T)
file=paste("data/",genomes[2],".df.",year,".blast.prot.tsv",sep="")
sacchC<-read.table(file, sep='\t', header=T)
file=paste("data/",genomes[3],".df.",year,".blast.prot.tsv",sep="")
#file=paste("data/","homo_sapiens.df.prot.tsv",sep="")
homoC<-read.table(file, sep='\t', header=T)


columns=NULL
cutoffs=NULL
cutoffs[1]=-1
cutoffs[2]=0
cutoffs[3]=100
cutoffs[4]=200
cutoffs[5]=500
cutoffs[6]=1000
cutoffs[7]=1000000000


PDBcutoffs=NULL
PDBcutoffs[1]=99999
PDBcutoffs[2]=1.e-100
PDBcutoffs[3]=1.e-40
PDBcutoffs[4]=1.e-10
PDBcutoffs[5]=1.e-4
PDBcutoffs[6]=1.e-2
PDBcutoffs[7]=-99999

                                        # ----- Pre processing ---

ecoliA$Pfam_Meff=ecoliA$Pfam_Meff100
ecoliA$Pfam_Meff[ecoliA$Pfam_Meff100>0]<-101
ecoliA$Pfam_Meff[ecoliA$Pfam_Meff200>0]<-201
ecoliA$Pfam_Meff[ecoliA$Pfam_Meff500>0]<-501
ecoliA$Pfam_Meff[ecoliA$Pfam_Meff1000>0]<-1001
ecoliA$PDB_ID=ecoliA$PDB_count
ecoliA$PDB_ID[ecoliA$PDB_count==0]<-""
ecoliA$PDB_ID[ecoliA$PDB_count>0]<-"PDB"

sacchA$Pfam_Meff=sacchA$Pfam_Meff100
sacchA$Pfam_Meff[sacchA$Pfam_Meff100>0]<-101
sacchA$Pfam_Meff[sacchA$Pfam_Meff200>0]<-201
sacchA$Pfam_Meff[sacchA$Pfam_Meff500>0]<-501
sacchA$Pfam_Meff[sacchA$Pfam_Meff1000>0]<-1001
sacchA$PDB_ID=sacchA$PDB_count
sacchA$PDB_ID[sacchA$PDB_count==0]<-""
sacchA$PDB_ID[sacchA$PDB_count>0]<-"PDB"

homoA$Pfam_Meff=homoA$Pfam_Meff100
homoA$Pfam_Meff[homoA$Pfam_Meff100>0]<-101
homoA$Pfam_Meff[homoA$Pfam_Meff200>0]<-201
homoA$Pfam_Meff[homoA$Pfam_Meff500>0]<-501
homoA$Pfam_Meff[homoA$Pfam_Meff1000>0]<-1001
homoA$PDB_ID=homoA$PDB_count
homoA$PDB_ID[homoA$PDB_count==0]<-""
homoA$PDB_ID[homoA$PDB_count>0]<-"PDB"

ecoliB$Pfam_Meff=ecoliB$Pfam_Meff100
ecoliB$Pfam_Meff[ecoliB$Pfam_Meff100>0]<-101
ecoliB$Pfam_Meff[ecoliB$Pfam_Meff200>0]<-201
ecoliB$Pfam_Meff[ecoliB$Pfam_Meff500>0]<-501
ecoliB$Pfam_Meff[ecoliB$Pfam_Meff1000>0]<-1001
ecoliB$PDB_ID=ecoliB$PDB_count
ecoliB$PDB_ID[ecoliB$PDB_count==0]<-""
ecoliB$PDB_ID[ecoliB$PDB_count>0]<-"PDB"

sacchB$Pfam_Meff=sacchB$Pfam_Meff100
sacchB$Pfam_Meff[sacchB$Pfam_Meff100>0]<-101
sacchB$Pfam_Meff[sacchB$Pfam_Meff200>0]<-201
sacchB$Pfam_Meff[sacchB$Pfam_Meff500>0]<-501
sacchB$Pfam_Meff[sacchB$Pfam_Meff1000>0]<-1001
sacchB$PDB_ID=sacchB$PDB_count
sacchB$PDB_ID[sacchB$PDB_count==0]<-""
sacchB$PDB_ID[sacchB$PDB_count>0]<-"PDB"

homoB$Pfam_Meff=homoB$Pfam_Meff100
homoB$Pfam_Meff[homoB$Pfam_Meff100>0]<-101
homoB$Pfam_Meff[homoB$Pfam_Meff200>0]<-201
homoB$Pfam_Meff[homoB$Pfam_Meff500>0]<-501
homoB$Pfam_Meff[homoB$Pfam_Meff1000>0]<-1001
homoB$PDB_ID=homoB$PDB_count
homoB$PDB_ID[homoB$PDB_count==0]<-""
homoB$PDB_ID[homoB$PDB_count>0]<-"PDB"

ecoliC$Pfam_Meff=ecoliC$Pfam_Meff100
ecoliC$Pfam_Meff[ecoliC$Pfam_Meff100>0]<-101
ecoliC$Pfam_Meff[ecoliC$Pfam_Meff200>0]<-201
ecoliC$Pfam_Meff[ecoliC$Pfam_Meff500>0]<-501
ecoliC$Pfam_Meff[ecoliC$Pfam_Meff1000>0]<-1001
ecoliC$PDB_ID=ecoliC$PDB_count
ecoliC$PDB_ID[ecoliC$PDB_count==0]<-""
ecoliC$PDB_ID[ecoliC$PDB_count>0]<-"PDB"

sacchC$Pfam_Meff=sacchC$Pfam_Meff100
sacchC$Pfam_Meff[sacchC$Pfam_Meff100>0]<-101
sacchC$Pfam_Meff[sacchC$Pfam_Meff200>0]<-201
sacchC$Pfam_Meff[sacchC$Pfam_Meff500>0]<-501
sacchC$Pfam_Meff[sacchC$Pfam_Meff1000>0]<-1001
sacchC$PDB_ID=sacchC$PDB_count
sacchC$PDB_ID[sacchC$PDB_count==0]<-""
sacchC$PDB_ID[sacchC$PDB_count>0]<-"PDB"

homoC$Pfam_Meff=homoC$Pfam_Meff100
homoC$Pfam_Meff[homoC$Pfam_Meff100>0]<-101
homoC$Pfam_Meff[homoC$Pfam_Meff200>0]<-201
homoC$Pfam_Meff[homoC$Pfam_Meff500>0]<-501
homoC$Pfam_Meff[homoC$Pfam_Meff1000>0]<-1001
homoC$PDB_ID=homoC$PDB_count
homoC$PDB_ID[homoC$PDB_count==0]<-""
homoC$PDB_ID[homoC$PDB_count>0]<-"PDB"



                                        #-----------------------  E. Coli  --------------------------------
ecoliAAll=NULL
ecoliATM=NULL
ecoliADiso=NULL
ecoliAPDB=NULL
ecoliASeg=NULL
ecoliANoTM=NULL
ecoliANoDiso=NULL
ecoliANoPDB=NULL
ecoliANoSeg=NULL
ecoliAPDBAll=NULL
ecoliAPDBTM=NULL
ecoliAPDBDiso=NULL
ecoliAPDBPDB=NULL
ecoliAPDBSeg=NULL
ecoliAPDBNoTM=NULL
ecoliAPDBNoDiso=NULL
ecoliAPDBNoSeg=NULL
ecoliAPDBTMNoDiso=NULL
ecoliAPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)


#for (i in loop){
#    ecoliAPDB[i]=length(which(ecoliA$PDB_E-value < cutoffs[i] & ecoliA$PDB_E-value >= PDBcutoffs[i+1] & ecoliA$PDB_ID != ""))
#}

                                        # Ja
for (i in loop){
    ecoliAAll[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1]))
    ecoliATM[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$TM > 0))
    ecoliANoTM[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$TM == 0))
    ecoliADiso[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$Disorder>=0.5))
    ecoliANoDiso[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$Disorder<0.5))
    ecoliAPDB[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$PDB_ID != ""))
    ecoliANoPDB[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$PDB_ID == ""))
    ecoliASeg[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$Seg_low > 0))
    ecoliANoSeg[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$Seg_low == 0))
    labels[i]=cutoffs[i]
    ecoliAPDBAll[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$PDB_ID == "") )
    ecoliAPDBTM[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$TM > 0 & ecoliA$PDB_ID == ""))
    ecoliAPDBNoTM[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$TM == 0  & ecoliA$PDB_ID == "" & ecoliA$Disorder<0.5 ))
    ecoliAPDBTMNoDiso[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$TM > 0 & ecoliA$PDB_ID == ""& ecoliA$Disorder<0.5 ))
    ecoliAPDBNoTMNoDiso[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$TM == 0  & ecoliA$PDB_ID == ""))
    ecoliAPDBDiso[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$Disorder>=0.5  & ecoliA$PDB_ID == ""))
    ecoliAPDBNoDiso[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$Disorder<0.5  & ecoliA$PDB_ID == ""))
    ecoliAPDBSeg[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$Seg_low > 0 & ecoliA$PDB_ID == ""))
    ecoliAPDBNoSeg[i]=length(which(ecoliA$Pfam_Meff > cutoffs[i] & ecoliA$Pfam_Meff <= cutoffs[i+1] & ecoliA$Seg_low == 0  & ecoliA$PDB_ID == ""))
}

ecoliAPDBAll[7]=length(which(ecoliA$PDB_ID != "") )
ecoliADisoAll=c(ecoliAPDBNoDiso,ecoliAPDBDiso,ecoliAPDBAll[7])
ecoliATMAll=c(ecoliAPDBNoTMNoDiso[1:6],ecoliAPDBTMNoDiso[1:6],ecoliAPDBDiso[1:6],ecoliAPDB[1:6])
ecoliAPDBTM[7]=length(which(ecoliA$PDB_ID != "" & ecoliA$TM > 0) )
ecoliAPDBNoTM[7]=length(which(ecoliA$PDB_ID != "" & ecoliA$TM == 0 ) )
ecoliAPDBDiso[7]=length(which(ecoliA$PDB_ID != "" & ecoliA$Disorder>=0.5 ) )
ecoliAPDBNoDiso[7]=length(which(ecoliA$PDB_ID != ""  & ecoliA$Disorder<0.5 ) )
ecoliAPDBSeg[7]=length(which(ecoliA$PDB_ID != "" & ecoliA$Seg_low > 0) )
ecoliAPDBNoSeg[7]=length(which(ecoliA$PDB_ID != ""  & ecoliA$Seg_low == 0) )


ecoliAAll[ecoliAAll==0]<-0.001
ecoliAPDBAll[ecoliAPDBAll==0]<-0.001
ecoliADisoAll[ecoliADisoAll==0]<-0.001
ecoliATMAll[ecoliATMAll==0]<-0.001

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


pct <- round(ecoliAAll/sum(ecoliAAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="ecoliA-prot"
outfile=paste("figures/",genome,".",year,"-PDB-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(ecoliAPDB)/(sum(ecoliANoPDB)+sum(ecoliAPDB)))**2
                                        #    pie(ecoliAAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliAPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliANoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliAAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliATM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliANoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliAAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliADiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliANoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="ecoliA-all-prot"
pct <- round(ecoliAPDBAll/sum(ecoliAAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliAPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliAPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliAPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliAPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliAPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliAPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()


                                        #-----------------------  E. Coli  --------------------------------
ecoliBAll=NULL
ecoliBTM=NULL
ecoliBDiso=NULL
ecoliBPDB=NULL
ecoliBSeg=NULL
ecoliBNoTM=NULL
ecoliBNoDiso=NULL
ecoliBNoPDB=NULL
ecoliBNoSeg=NULL
ecoliBPDBAll=NULL
ecoliBPDBTM=NULL
ecoliBPDBDiso=NULL
ecoliBPDBPDB=NULL
ecoliBPDBSeg=NULL
ecoliBPDBNoTM=NULL
ecoliBPDBNoDiso=NULL
ecoliBPDBNoSeg=NULL
ecoliBPDBTMNoDiso=NULL
ecoliBPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)




                                        # Ja
for (i in loop){
    ecoliBAll[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1]))
    ecoliBTM[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$TM > 0))
    ecoliBNoTM[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$TM == 0))
    ecoliBDiso[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$Disorder>=0.5))
    ecoliBNoDiso[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$Disorder<0.5))
    ecoliBPDB[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$PDB_ID != ""))
    ecoliBNoPDB[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$PDB_ID == ""))
    ecoliBSeg[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$Seg_low > 0))
    ecoliBNoSeg[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$Seg_low == 0))
    labels[i]=cutoffs[i]
    ecoliBPDBAll[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$PDB_ID == "") )
    ecoliBPDBTM[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$TM > 0 & ecoliB$PDB_ID == ""))
    ecoliBPDBNoTM[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$TM == 0  & ecoliB$PDB_ID == "" & ecoliB$Disorder<0.5 ))
    ecoliBPDBTMNoDiso[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$TM > 0 & ecoliB$PDB_ID == ""& ecoliB$Disorder<0.5 ))
    ecoliBPDBNoTMNoDiso[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$TM == 0  & ecoliB$PDB_ID == ""))
    ecoliBPDBDiso[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$Disorder>=0.5  & ecoliB$PDB_ID == ""))
    ecoliBPDBNoDiso[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$Disorder<0.5  & ecoliB$PDB_ID == ""))
    ecoliBPDBSeg[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$Seg_low > 0 & ecoliB$PDB_ID == ""))
    ecoliBPDBNoSeg[i]=length(which(ecoliB$Pfam_Meff > cutoffs[i] & ecoliB$Pfam_Meff <= cutoffs[i+1] & ecoliB$Seg_low == 0  & ecoliB$PDB_ID == ""))
}

ecoliBPDBAll[7]=length(which(ecoliB$PDB_ID != "") )
ecoliBDisoAll=c(ecoliBPDBNoDiso,ecoliBPDBDiso,ecoliBPDBAll[7])
ecoliBTMAll=c(ecoliBPDBNoTMNoDiso[1:6],ecoliBPDBTMNoDiso[1:6],ecoliBPDBDiso[1:6],ecoliBPDB[1:6])
ecoliBPDBTM[7]=length(which(ecoliB$PDB_ID != "" & ecoliB$TM > 0) )
ecoliBPDBNoTM[7]=length(which(ecoliB$PDB_ID != "" & ecoliB$TM == 0 ) )
ecoliBPDBDiso[7]=length(which(ecoliB$PDB_ID != "" & ecoliB$Disorder>=0.5 ) )
ecoliBPDBNoDiso[7]=length(which(ecoliB$PDB_ID != ""  & ecoliB$Disorder<0.5 ) )
ecoliBPDBSeg[7]=length(which(ecoliB$PDB_ID != "" & ecoliB$Seg_low > 0) )
ecoliBPDBNoSeg[7]=length(which(ecoliB$PDB_ID != ""  & ecoliB$Seg_low == 0) )


ecoliBAll[ecoliBAll==0]<-0.001
ecoliBPDBAll[ecoliBPDBAll==0]<-0.001
ecoliBDisoAll[ecoliBDisoAll==0]<-0.001
ecoliBTMAll[ecoliBTMAll==0]<-0.001

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


pct <- round(ecoliBAll/sum(ecoliBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="ecoliB-prot"
outfile=paste("figures/",genome,".",year,"-PDB-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(ecoliBPDB)/(sum(ecoliBNoPDB)+sum(ecoliBPDB)))**2
                                        #    pie(ecoliBAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliBPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliBNoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="ecoliB-all-prot"
pct <- round(ecoliBPDBAll/sum(ecoliBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliBPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliBPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliBPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliBPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliBPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliBPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()


                                        #-----------------------  E. Coli  --------------------------------
ecoliCAll=NULL
ecoliCTM=NULL
ecoliCDiso=NULL
ecoliCPDB=NULL
ecoliCSeg=NULL
ecoliCNoTM=NULL
ecoliCNoDiso=NULL
ecoliCNoPDB=NULL
ecoliCNoSeg=NULL
ecoliCPDBAll=NULL
ecoliCPDBTM=NULL
ecoliCPDBDiso=NULL
ecoliCPDBPDB=NULL
ecoliCPDBSeg=NULL
ecoliCPDBNoTM=NULL
ecoliCPDBNoDiso=NULL
ecoliCPDBNoSeg=NULL
ecoliCPDBTMNoDiso=NULL
ecoliCPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)




                                        # Ja
for (i in loop){
    ecoliCAll[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1]))
    ecoliCTM[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$TM > 0))
    ecoliCNoTM[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$TM == 0))
    ecoliCDiso[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$Disorder>=0.5))
    ecoliCNoDiso[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$Disorder<0.5))
    ecoliCPDB[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$PDB_ID != ""))
    ecoliCNoPDB[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$PDB_ID == ""))
    ecoliCSeg[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$Seg_low > 0))
    ecoliCNoSeg[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$Seg_low == 0))
    labels[i]=cutoffs[i]
    ecoliCPDBAll[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$PDB_ID == "") )
    ecoliCPDBTM[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$TM > 0 & ecoliC$PDB_ID == ""))
    ecoliCPDBNoTM[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$TM == 0  & ecoliC$PDB_ID == "" & ecoliC$Disorder<0.5 ))
    ecoliCPDBTMNoDiso[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$TM > 0 & ecoliC$PDB_ID == ""& ecoliC$Disorder<0.5 ))
    ecoliCPDBNoTMNoDiso[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$TM == 0  & ecoliC$PDB_ID == ""))
    ecoliCPDBDiso[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$Disorder>=0.5  & ecoliC$PDB_ID == ""))
    ecoliCPDBNoDiso[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$Disorder<0.5  & ecoliC$PDB_ID == ""))
    ecoliCPDBSeg[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$Seg_low > 0 & ecoliC$PDB_ID == ""))
    ecoliCPDBNoSeg[i]=length(which(ecoliC$Pfam_Meff > cutoffs[i] & ecoliC$Pfam_Meff <= cutoffs[i+1] & ecoliC$Seg_low == 0  & ecoliC$PDB_ID == ""))
}

ecoliCPDBAll[7]=length(which(ecoliC$PDB_ID != "") )
ecoliCDisoAll=c(ecoliCPDBNoDiso,ecoliCPDBDiso,ecoliCPDBAll[7])
ecoliCTMAll=c(ecoliCPDBNoTMNoDiso[1:6],ecoliCPDBTMNoDiso[1:6],ecoliCPDBDiso[1:6],ecoliCPDB[1:6])
ecoliCPDBTM[7]=length(which(ecoliC$PDB_ID != "" & ecoliC$TM > 0) )
ecoliCPDBNoTM[7]=length(which(ecoliC$PDB_ID != "" & ecoliC$TM == 0 ) )
ecoliCPDBDiso[7]=length(which(ecoliC$PDB_ID != "" & ecoliC$Disorder>=0.5 ) )
ecoliCPDBNoDiso[7]=length(which(ecoliC$PDB_ID != ""  & ecoliC$Disorder<0.5 ) )
ecoliCPDBSeg[7]=length(which(ecoliC$PDB_ID != "" & ecoliC$Seg_low > 0) )
ecoliCPDBNoSeg[7]=length(which(ecoliC$PDB_ID != ""  & ecoliC$Seg_low == 0) )


ecoliCAll[ecoliCAll==0]<-0.001
ecoliCPDBAll[ecoliCPDBAll==0]<-0.001
ecoliCDisoAll[ecoliCDisoAll==0]<-0.001
ecoliCTMAll[ecoliCTMAll==0]<-0.001

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


pct <- round(ecoliCAll/sum(ecoliCAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="ecoliC-prot"
outfile=paste("figures/",genome,".",year,"-PDB-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(ecoliCPDB)/(sum(ecoliCNoPDB)+sum(ecoliCPDB)))**2
                                        #    pie(ecoliCAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliCPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliCNoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliCAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliCTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliCNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliCAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliCDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliCNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="ecoliC-all-prot"
pct <- round(ecoliCPDBAll/sum(ecoliCAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliCPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliCPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliCPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(ecoliCPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,ecoliCPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliCPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()


                                        #-----------------------  Sacharomyces  --------------------------------
sacchAAll=NULL
sacchATM=NULL
sacchADiso=NULL
sacchAPDB=NULL
sacchASeg=NULL
sacchANoTM=NULL
sacchANoDiso=NULL
sacchANoPDB=NULL
sacchANoSeg=NULL
sacchAPDBAll=NULL
sacchAPDBTM=NULL
sacchAPDBDiso=NULL
sacchAPDBPDB=NULL
sacchAPDBSeg=NULL
sacchAPDBNoTM=NULL
sacchAPDBNoDiso=NULL
sacchAPDBNoSeg=NULL
sacchAPDBTMNoDiso=NULL
sacchAPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)




                                        # Ja
for (i in loop){
    sacchAAll[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1]))
    sacchATM[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$TM > 0))
    sacchANoTM[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$TM == 0))
    sacchADiso[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$Disorder>=0.5))
    sacchANoDiso[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$Disorder<0.5))
    sacchAPDB[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$PDB_ID != ""))
    sacchANoPDB[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$PDB_ID == ""))
    sacchASeg[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$Seg_low > 0))
    sacchANoSeg[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$Seg_low == 0))
    labels[i]=cutoffs[i]
    sacchAPDBAll[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$PDB_ID == "") )
    sacchAPDBTM[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$TM > 0 & sacchA$PDB_ID == ""))
    sacchAPDBNoTM[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$TM == 0  & sacchA$PDB_ID == ""))
    sacchAPDBTMNoDiso[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$TM > 0 & sacchA$PDB_ID == ""& sacchA$Disorder<0.5 ))
    sacchAPDBNoTMNoDiso[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$TM == 0  & sacchA$PDB_ID == "" & sacchA$Disorder<0.5 ))
    sacchAPDBDiso[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$Disorder>=0.5  & sacchA$PDB_ID == ""))
    sacchAPDBNoDiso[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$Disorder<0.5  & sacchA$PDB_ID == ""))
    sacchAPDBSeg[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$Seg_low > 0 & sacchA$PDB_ID == ""))
    sacchAPDBNoSeg[i]=length(which(sacchA$Pfam_Meff > cutoffs[i] & sacchA$Pfam_Meff <= cutoffs[i+1] & sacchA$Seg_low == 0  & sacchA$PDB_ID == ""))
}
sacchAPDBAll[7]=length(which(sacchA$PDB_ID != "") )
sacchATMAll=c(sacchAPDBNoTMNoDiso[1:6],sacchAPDBTMNoDiso[1:6],sacchAPDBDiso[1:6],sacchAPDB[1:6])
sacchADisoAll=c(sacchAPDBNoDiso,sacchAPDBDiso,sacchAPDBAll[7])
sacchAPDBTM[7]=length(which(sacchA$PDB_ID != "" & sacchA$TM > 0) )
sacchAPDBNoTM[7]=length(which(sacchA$PDB_ID != "" & sacchA$TM == 0 ) )
sacchAPDBDiso[7]=length(which(sacchA$PDB_ID != "" & sacchA$Disorder>=0.5 ) )
sacchAPDBNoDiso[7]=length(which(sacchA$PDB_ID != ""  & sacchA$Disorder<0.5 ) )
sacchAPDBSeg[7]=length(which(sacchA$PDB_ID != "" & sacchA$Seg_low > 0) )
sacchAPDBNoSeg[7]=length(which(sacchA$PDB_ID != ""  & sacchA$Seg_low == 0) )


sacchAAll[sacchAAll==0]<-0.001
sacchAPDBAll[sacchAPDBAll==0]<-0.001
sacchADisoAll[sacchADisoAll==0]<-0.001
sacchATMAll[sacchATMAll==0]<-0.001

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


pct <- round(sacchAAll/sum(sacchAAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="sacchA-prot"
outfile=paste("figures/",genome,".",year,"-PDB.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(sacchAPDB)/(sum(sacchANoPDB)+sum(sacchAPDB)))**2
                                        #    pie(sacchAAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchAPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchANoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchAAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchATM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchANoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchAAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchADiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchANoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="sacchA-all-prot"
pct <- round(sacchAAll/sum(sacchAPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchAPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchAPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchAPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchAPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchAPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchAPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()



                                        #-----------------------  Sacharomyces  --------------------------------
sacchBAll=NULL
sacchBTM=NULL
sacchBDiso=NULL
sacchBPDB=NULL
sacchBSeg=NULL
sacchBNoTM=NULL
sacchBNoDiso=NULL
sacchBNoPDB=NULL
sacchBNoSeg=NULL
sacchBPDBAll=NULL
sacchBPDBTM=NULL
sacchBPDBDiso=NULL
sacchBPDBPDB=NULL
sacchBPDBSeg=NULL
sacchBPDBNoTM=NULL
sacchBPDBNoDiso=NULL
sacchBPDBNoSeg=NULL
sacchBPDBTMNoDiso=NULL
sacchBPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)




                                        # Ja
for (i in loop){
    sacchBAll[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1]))
    sacchBTM[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$TM > 0))
    sacchBNoTM[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$TM == 0))
    sacchBDiso[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$Disorder>=0.5))
    sacchBNoDiso[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$Disorder<0.5))
    sacchBPDB[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$PDB_ID != ""))
    sacchBNoPDB[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$PDB_ID == ""))
    sacchBSeg[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$Seg_low > 0))
    sacchBNoSeg[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$Seg_low == 0))
    labels[i]=cutoffs[i]
    sacchBPDBAll[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$PDB_ID == "") )
    sacchBPDBTM[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$TM > 0 & sacchB$PDB_ID == ""))
    sacchBPDBNoTM[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$TM == 0  & sacchB$PDB_ID == ""))
    sacchBPDBTMNoDiso[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$TM > 0 & sacchB$PDB_ID == ""& sacchB$Disorder<0.5 ))
    sacchBPDBNoTMNoDiso[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$TM == 0  & sacchB$PDB_ID == "" & sacchB$Disorder<0.5 ))
    sacchBPDBDiso[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$Disorder>=0.5  & sacchB$PDB_ID == ""))
    sacchBPDBNoDiso[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$Disorder<0.5  & sacchB$PDB_ID == ""))
    sacchBPDBSeg[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$Seg_low > 0 & sacchB$PDB_ID == ""))
    sacchBPDBNoSeg[i]=length(which(sacchB$Pfam_Meff > cutoffs[i] & sacchB$Pfam_Meff <= cutoffs[i+1] & sacchB$Seg_low == 0  & sacchB$PDB_ID == ""))
}
sacchBPDBAll[7]=length(which(sacchB$PDB_ID != "") )
sacchBTMAll=c(sacchBPDBNoTMNoDiso[1:6],sacchBPDBTMNoDiso[1:6],sacchBPDBDiso[1:6],sacchBPDB[1:6])
sacchBDisoAll=c(sacchBPDBNoDiso,sacchBPDBDiso,sacchBPDBAll[7])
sacchBPDBTM[7]=length(which(sacchB$PDB_ID != "" & sacchB$TM > 0) )
sacchBPDBNoTM[7]=length(which(sacchB$PDB_ID != "" & sacchB$TM == 0 ) )
sacchBPDBDiso[7]=length(which(sacchB$PDB_ID != "" & sacchB$Disorder>=0.5 ) )
sacchBPDBNoDiso[7]=length(which(sacchB$PDB_ID != ""  & sacchB$Disorder<0.5 ) )
sacchBPDBSeg[7]=length(which(sacchB$PDB_ID != "" & sacchB$Seg_low > 0) )
sacchBPDBNoSeg[7]=length(which(sacchB$PDB_ID != ""  & sacchB$Seg_low == 0) )


sacchBAll[sacchBAll==0]<-0.001
sacchBPDBAll[sacchBPDBAll==0]<-0.001
sacchBDisoAll[sacchBDisoAll==0]<-0.001
sacchBTMAll[sacchBTMAll==0]<-0.001

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


pct <- round(sacchBAll/sum(sacchBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="sacchB-prot"
outfile=paste("figures/",genome,".",year,"-PDB.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(sacchBPDB)/(sum(sacchBNoPDB)+sum(sacchBPDB)))**2
                                        #    pie(sacchBAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchBPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchBNoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="sacchB-all-prot"
pct <- round(sacchBAll/sum(sacchBPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchBPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchBPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchBPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchBPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchBPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchBPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()



                                        #-----------------------  Sacharomyces  --------------------------------
sacchCAll=NULL
sacchCTM=NULL
sacchCDiso=NULL
sacchCPDB=NULL
sacchCSeg=NULL
sacchCNoTM=NULL
sacchCNoDiso=NULL
sacchCNoPDB=NULL
sacchCNoSeg=NULL
sacchCPDBAll=NULL
sacchCPDBTM=NULL
sacchCPDBDiso=NULL
sacchCPDBPDB=NULL
sacchCPDBSeg=NULL
sacchCPDBNoTM=NULL
sacchCPDBNoDiso=NULL
sacchCPDBNoSeg=NULL
sacchCPDBTMNoDiso=NULL
sacchCPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)




                                        # Ja
for (i in loop){
    sacchCAll[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1]))
    sacchCTM[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$TM > 0))
    sacchCNoTM[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$TM == 0))
    sacchCDiso[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$Disorder>=0.5))
    sacchCNoDiso[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$Disorder<0.5))
    sacchCPDB[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$PDB_ID != ""))
    sacchCNoPDB[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$PDB_ID == ""))
    sacchCSeg[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$Seg_low > 0))
    sacchCNoSeg[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$Seg_low == 0))
    labels[i]=cutoffs[i]
    sacchCPDBAll[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$PDB_ID == "") )
    sacchCPDBTM[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$TM > 0 & sacchC$PDB_ID == ""))
    sacchCPDBNoTM[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$TM == 0  & sacchC$PDB_ID == ""))
    sacchCPDBTMNoDiso[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$TM > 0 & sacchC$PDB_ID == ""& sacchC$Disorder<0.5 ))
    sacchCPDBNoTMNoDiso[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$TM == 0  & sacchC$PDB_ID == "" & sacchC$Disorder<0.5 ))
    sacchCPDBDiso[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$Disorder>=0.5  & sacchC$PDB_ID == ""))
    sacchCPDBNoDiso[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$Disorder<0.5  & sacchC$PDB_ID == ""))
    sacchCPDBSeg[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$Seg_low > 0 & sacchC$PDB_ID == ""))
    sacchCPDBNoSeg[i]=length(which(sacchC$Pfam_Meff > cutoffs[i] & sacchC$Pfam_Meff <= cutoffs[i+1] & sacchC$Seg_low == 0  & sacchC$PDB_ID == ""))
}
sacchCPDBAll[7]=length(which(sacchC$PDB_ID != "") )
sacchCTMAll=c(sacchCPDBNoTMNoDiso[1:6],sacchCPDBTMNoDiso[1:6],sacchCPDBDiso[1:6],sacchCPDB[1:6])
sacchCDisoAll=c(sacchCPDBNoDiso,sacchCPDBDiso,sacchCPDBAll[7])
sacchCPDBTM[7]=length(which(sacchC$PDB_ID != "" & sacchC$TM > 0) )
sacchCPDBNoTM[7]=length(which(sacchC$PDB_ID != "" & sacchC$TM == 0 ) )
sacchCPDBDiso[7]=length(which(sacchC$PDB_ID != "" & sacchC$Disorder>=0.5 ) )
sacchCPDBNoDiso[7]=length(which(sacchC$PDB_ID != ""  & sacchC$Disorder<0.5 ) )
sacchCPDBSeg[7]=length(which(sacchC$PDB_ID != "" & sacchC$Seg_low > 0) )
sacchCPDBNoSeg[7]=length(which(sacchC$PDB_ID != ""  & sacchC$Seg_low == 0) )


sacchCAll[sacchCAll==0]<-0.001
sacchCPDBAll[sacchCPDBAll==0]<-0.001
sacchCDisoAll[sacchCDisoAll==0]<-0.001
sacchCTMAll[sacchCTMAll==0]<-0.001

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


pct <- round(sacchCAll/sum(sacchCAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="sacchC-prot"
outfile=paste("figures/",genome,".",year,"-PDB.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(sacchCPDB)/(sum(sacchCNoPDB)+sum(sacchCPDB)))**2
                                        #    pie(sacchCAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchCPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchCNoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchCAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchCTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchCNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchCAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchCDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchCNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="sacchC-all-prot"
pct <- round(sacchCAll/sum(sacchCPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchCPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchCPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchCPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(sacchCPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,sacchCPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,sacchCPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()



                                        #-----------------------  Homo Sapiens  --------------------------------
homoAAll=NULL
homoATM=NULL
homoADiso=NULL
homoAPDB=NULL
homoASeg=NULL
homoANoTM=NULL
homoANoDiso=NULL
homoANoPDB=NULL
homoANoSeg=NULL
homoAPDBAll=NULL
homoAPDBTM=NULL
homoAPDBDiso=NULL
homoAPDBPDB=NULL
homoAPDBSeg=NULL
homoAPDBNoTM=NULL
homoAPDBNoDiso=NULL
homoAPDBNoSeg=NULL
homoAPDBTMNoDiso=NULL
homoAPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)





                                        # Ja
for (i in loop){
    homoAAll[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1]))
    homoATM[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$TM > 0))
    homoANoTM[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$TM == 0))
    homoADiso[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$Disorder>=0.5))
    homoANoDiso[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$Disorder<0.5))
    homoAPDB[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$PDB_ID != ""))
    homoANoPDB[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$PDB_ID == ""))
    homoASeg[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$Seg_low > 0))
    homoANoSeg[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$Seg_low == 0))
    labels[i]=cutoffs[i]
    homoAPDBAll[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$PDB_ID == "") )
    homoAPDBTM[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$TM > 0 & homoA$PDB_ID == ""))
    homoAPDBNoTM[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$TM == 0  & homoA$PDB_ID == ""))
    homoAPDBNoTMNoDiso[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$TM == 0  & homoA$PDB_ID == "" & homoA$Disorder<0.5 ))
    homoAPDBTMNoDiso[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$TM > 0 & homoA$PDB_ID == ""& homoA$Disorder<0.5 ))
    homoAPDBDiso[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$Disorder>=0.5  & homoA$PDB_ID == ""))
    homoAPDBNoDiso[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$Disorder<0.5  & homoA$PDB_ID == ""))
    homoAPDBSeg[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$Seg_low > 0 & homoA$PDB_ID == ""))
    homoAPDBNoSeg[i]=length(which(homoA$Pfam_Meff > cutoffs[i] & homoA$Pfam_Meff <= cutoffs[i+1] & homoA$Seg_low == 0  & homoA$PDB_ID == ""))
}
homoAPDBAll[7]=length(which(homoA$PDB_ID != "") )
homoADisoAll=c(homoAPDBNoDiso,homoAPDBDiso,homoAPDBAll[7])
homoATMAll=c(homoAPDBNoTMNoDiso[1:6],homoAPDBTMNoDiso[1:6],homoAPDBDiso[1:6],homoAPDB[1:6])
homoAPDBTM[7]=length(which(homoA$PDB_ID != "" & homoA$TM > 0) )
homoAPDBNoTM[7]=length(which(homoA$PDB_ID != "" & homoA$TM == 0 ) )
homoAPDBDiso[7]=length(which(homoA$PDB_ID != "" & homoA$Disorder>=0.5 ) )
homoAPDBNoDiso[7]=length(which(homoA$PDB_ID != ""  & homoA$Disorder<0.5 ) )
homoAPDBSeg[7]=length(which(homoA$PDB_ID != "" & homoA$Seg_low > 0) )
homoAPDBNoSeg[7]=length(which(homoA$PDB_ID != ""  & homoA$Seg_low == 0) )

homoAAll[homoAAll==0]<-0.001
homoAPDBAll[homoAPDBAll==0]<-0.001
homoADisoAll[homoADisoAll==0]<-0.001
homoATMAll[homoATMAll==0]<-0.001

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


pct <- round(homoAAll/sum(homoAAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="homoA-prot"
outfile=paste("figures/",genome,".",year,"-PDB.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(homoAPDB)/(sum(homoANoPDB)+sum(homoAPDB)))**2
                                        #    pie(homoAAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoAPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoANoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoAAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoATM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoANoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoAAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoADiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoANoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="homoA-all-prot"
pct <- round(homoAAll/sum(homoAPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoAPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoAPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoAPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoAPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoAPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoAPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()


                                        #-----------------------  Homo Sapiens  --------------------------------
homoBAll=NULL
homoBTM=NULL
homoBDiso=NULL
homoBPDB=NULL
homoBSeg=NULL
homoBNoTM=NULL
homoBNoDiso=NULL
homoBNoPDB=NULL
homoBNoSeg=NULL
homoBPDBAll=NULL
homoBPDBTM=NULL
homoBPDBDiso=NULL
homoBPDBPDB=NULL
homoBPDBSeg=NULL
homoBPDBNoTM=NULL
homoBPDBNoDiso=NULL
homoBPDBNoSeg=NULL
homoBPDBTMNoDiso=NULL
homoBPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)





                                        # Ja
for (i in loop){
    homoBAll[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1]))
    homoBTM[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$TM > 0))
    homoBNoTM[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$TM == 0))
    homoBDiso[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$Disorder>=0.5))
    homoBNoDiso[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$Disorder<0.5))
    homoBPDB[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$PDB_ID != ""))
    homoBNoPDB[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$PDB_ID == ""))
    homoBSeg[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$Seg_low > 0))
    homoBNoSeg[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$Seg_low == 0))
    labels[i]=cutoffs[i]
    homoBPDBAll[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$PDB_ID == "") )
    homoBPDBTM[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$TM > 0 & homoB$PDB_ID == ""))
    homoBPDBNoTM[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$TM == 0  & homoB$PDB_ID == ""))
    homoBPDBNoTMNoDiso[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$TM == 0  & homoB$PDB_ID == "" & homoB$Disorder<0.5 ))
    homoBPDBTMNoDiso[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$TM > 0 & homoB$PDB_ID == ""& homoB$Disorder<0.5 ))
    homoBPDBDiso[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$Disorder>=0.5  & homoB$PDB_ID == ""))
    homoBPDBNoDiso[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$Disorder<0.5  & homoB$PDB_ID == ""))
    homoBPDBSeg[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$Seg_low > 0 & homoB$PDB_ID == ""))
    homoBPDBNoSeg[i]=length(which(homoB$Pfam_Meff > cutoffs[i] & homoB$Pfam_Meff <= cutoffs[i+1] & homoB$Seg_low == 0  & homoB$PDB_ID == ""))
}
homoBPDBAll[7]=length(which(homoB$PDB_ID != "") )
homoBDisoAll=c(homoBPDBNoDiso,homoBPDBDiso,homoBPDBAll[7])
homoBTMAll=c(homoBPDBNoTMNoDiso[1:6],homoBPDBTMNoDiso[1:6],homoBPDBDiso[1:6],homoBPDB[1:6])
homoBPDBTM[7]=length(which(homoB$PDB_ID != "" & homoB$TM > 0) )
homoBPDBNoTM[7]=length(which(homoB$PDB_ID != "" & homoB$TM == 0 ) )
homoBPDBDiso[7]=length(which(homoB$PDB_ID != "" & homoB$Disorder>=0.5 ) )
homoBPDBNoDiso[7]=length(which(homoB$PDB_ID != ""  & homoB$Disorder<0.5 ) )
homoBPDBSeg[7]=length(which(homoB$PDB_ID != "" & homoB$Seg_low > 0) )
homoBPDBNoSeg[7]=length(which(homoB$PDB_ID != ""  & homoB$Seg_low == 0) )

homoBAll[homoBAll==0]<-0.001
homoBPDBAll[homoBPDBAll==0]<-0.001
homoBDisoAll[homoBDisoAll==0]<-0.001
homoBTMAll[homoBTMAll==0]<-0.001

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


pct <- round(homoBAll/sum(homoBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="homoB-prot"
outfile=paste("figures/",genome,".",year,"-PDB.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(homoBPDB)/(sum(homoBNoPDB)+sum(homoBPDB)))**2
                                        #    pie(homoBAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoBPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoBNoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="homoB-all-prot"
pct <- round(homoBAll/sum(homoBPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoBPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoBPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoBPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoBPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoBPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoBPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()





                                        #-----------------------  Homo Sapiens  --------------------------------
homoCAll=NULL
homoCTM=NULL
homoCDiso=NULL
homoCPDB=NULL
homoCSeg=NULL
homoCNoTM=NULL
homoCNoDiso=NULL
homoCNoPDB=NULL
homoCNoSeg=NULL
homoCPDBAll=NULL
homoCPDBTM=NULL
homoCPDBDiso=NULL
homoCPDBPDB=NULL
homoCPDBSeg=NULL
homoCPDBNoTM=NULL
homoCPDBNoDiso=NULL
homoCPDBNoSeg=NULL
homoCPDBTMNoDiso=NULL
homoCPDBNoTMNoDiso=NULL
labels=NULL
loop=seq(1,6)





                                        # Ja
for (i in loop){
    homoCAll[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1]))
    homoCTM[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$TM > 0))
    homoCNoTM[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$TM == 0))
    homoCDiso[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$Disorder>=0.5))
    homoCNoDiso[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$Disorder<0.5))
    homoCPDB[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$PDB_ID != ""))
    homoCNoPDB[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$PDB_ID == ""))
    homoCSeg[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$Seg_low > 0))
    homoCNoSeg[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$Seg_low == 0))
    labels[i]=cutoffs[i]
    homoCPDBAll[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$PDB_ID == "") )
    homoCPDBTM[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$TM > 0 & homoC$PDB_ID == ""))
    homoCPDBNoTM[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$TM == 0  & homoC$PDB_ID == ""))
    homoCPDBNoTMNoDiso[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$TM == 0  & homoC$PDB_ID == "" & homoC$Disorder<0.5 ))
    homoCPDBTMNoDiso[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$TM > 0 & homoC$PDB_ID == ""& homoC$Disorder<0.5 ))
    homoCPDBDiso[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$Disorder>=0.5  & homoC$PDB_ID == ""))
    homoCPDBNoDiso[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$Disorder<0.5  & homoC$PDB_ID == ""))
    homoCPDBSeg[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$Seg_low > 0 & homoC$PDB_ID == ""))
    homoCPDBNoSeg[i]=length(which(homoC$Pfam_Meff > cutoffs[i] & homoC$Pfam_Meff <= cutoffs[i+1] & homoC$Seg_low == 0  & homoC$PDB_ID == ""))
}
homoCPDBAll[7]=length(which(homoC$PDB_ID != "") )
homoCDisoAll=c(homoCPDBNoDiso,homoCPDBDiso,homoCPDBAll[7])
homoCTMAll=c(homoCPDBNoTMNoDiso[1:6],homoCPDBTMNoDiso[1:6],homoCPDBDiso[1:6],homoCPDB[1:6])
homoCPDBTM[7]=length(which(homoC$PDB_ID != "" & homoC$TM > 0) )
homoCPDBNoTM[7]=length(which(homoC$PDB_ID != "" & homoC$TM == 0 ) )
homoCPDBDiso[7]=length(which(homoC$PDB_ID != "" & homoC$Disorder>=0.5 ) )
homoCPDBNoDiso[7]=length(which(homoC$PDB_ID != ""  & homoC$Disorder<0.5 ) )
homoCPDBSeg[7]=length(which(homoC$PDB_ID != "" & homoC$Seg_low > 0) )
homoCPDBNoSeg[7]=length(which(homoC$PDB_ID != ""  & homoC$Seg_low == 0) )

homoCAll[homoCAll==0]<-0.001
homoCPDBAll[homoCPDBAll==0]<-0.001
homoCDisoAll[homoCDisoAll==0]<-0.001
homoCTMAll[homoCTMAll==0]<-0.001

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


pct <- round(homoCAll/sum(homoCAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")
iniR=1
genome="homoC-prot"
outfile=paste("figures/",genome,".",year,"-PDB.png",sep="")
                                        #png(outfile,width=1280,height=1280)
fraction=(sum(homoCPDB)/(sum(homoCNoPDB)+sum(homoCPDB)))**2
                                        #    pie(homoCAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoCPDB, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoCNoPDB, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"PDB",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoCAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoCTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoCNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso1-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoCAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoCDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoCNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()

labels[7]="PDB"
colors[7]="green"

genome="homoC-all-prot"
pct <- round(homoCAll/sum(homoCPDBAll)*100,digits=1)
names <- paste(labels,pct)
names <- paste(names,"%",sep="")

outfile=paste("figures/",genome,".",year,"-TM-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoCPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoCPDBTM, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoCPDBNoTM, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"TM",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoTM",cex=2.0,border=NA)
                                        #dev.off()

outfile=paste("figures/",genome,".",year,"-Diso2-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoCPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
                                        #    floating.pie(0,0,homoCPDBDiso, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,homoCPDBNoDiso, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
##legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
                                        #legend(-1,.0,"all",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Diso",cex=2.0,border=NA)
                                        #legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
                                        #dev.off()


                                        #------------------------------------------------------------------------
                                        #Summary of all
genome="All-prot"
outfile=paste("figures/",genome,".",year,"-summary-pie.png",sep="")
                                        #png(outfile,width=1280,height=1280)

                                        #    pie(homoPDBAll, labels=labels,col=colors,sub="Genome Classification",main=year,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
                                        #    floating.pie(0,0,sacchPDBAll, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliPDBAll, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"Homo",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Sacc",cex=2.0,border=NA)
                                        #legend(-.4,.0,"EColi",cex=2.0,border=NA)
                                        #dev.off()


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
colors[12]="deeppink4"
labels[13]="PDB"
colors[13]="green"

                                        #ecoliDisoAll=c(ecoliPDBNoDiso,ecoliPDBDiso,ecoliPDBAll[7])
                                        #sacchDisoAll=c(sacchPDBNoDiso,sacchPDBDiso,sacchPDBAll[7])
                                        #homoDisoAll=c(homoPDBNoDiso,homoPDBDiso,homoPDBAll[7])

genome="All-prot"
outfile=paste("figures/",genome,".",year,"-Detailed1-pie.png",sep="")

pct <- round(ecoliDisoAll/sum(ecoliDisoAll)*100,digits=1)
Names <- paste(labels,pct)
names <- paste(names,"%",sep="")
                                        #png(outfile,width=1280,height=1280)

                                        #    pie(homoDisoAll, labels=labels,col=colors,sub="Detailed Genome Classification",main=year,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
                                        #    floating.pie(0,0,sacchDisoAll, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliDisoAll, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"Homo",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Sacch",cex=2.0,border=NA)
                                        #legend(-.4,.0,"EColi",cex=2.0,border=NA)
                                        #dev.off()

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
Labels[2]=" "
Labels[3]=" "
Labels[4]=" "
Labels[5]=" "
Labels[6]=" "
Labels[7]="TM"
Labels[8]=" "
Labels[9]=" "
Labels[10]=" "
Labels[11]=" "
Labels[12]=" "
Labels[13]="Disorder"
Labels[14]=" "
Labels[15]=" "
Labels[16]=" "
Labels[17]=" "
Labels[18]=" "
Labels[19]="PDB"
Labels[21]=" "
Labels[22]=" "
Labels[23]=" "
Labels[24]=" "


                                        #png(outfile,width=1280,height=1280)
                                        #    pie(homoTMAll, labels=labels,col=colors,sub="Genome Classification",main=year,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
                                        #    floating.pie(0,0,sacchTMAll, col=colors,main='',radius=0.7)
                                        #    floating.pie(0,0,ecoliTMAll, col=colors,main='',radius=0.4)
                                        #    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(-1,.0,"Homo",cex=2.0,border=NA)
                                        #legend(-.7,.0,"Sacch",cex=2.0,border=NA)
                                        #legend(-.4,.0,"EColi",cex=2.0,border=NA)
                                        #dev.off()

proteinsA=NULL
proteinsA$ecoli=NULL
proteinsA$count=NULL
proteinsA$count[1]=1000
proteinsA$count[2]=500
proteinsA$count[3]=200
proteinsA$count[4]=100
proteinsA$ecoliTM[1]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff1000>0 & ecoliA$TM==1))
proteinsA$ecoliTM[2]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff500>0 & ecoliA$TM==1))
proteinsA$ecoliTM[3]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff200>0 & ecoliA$TM==1))
proteinsA$ecoliTM[4]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff100>0 & ecoliA$TM==1))
proteinsA$ecoliNoTM[1]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff1000>0 & ecoliA$TM==0))
proteinsA$ecoliNoTM[2]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff500>0 & ecoliA$TM==0))
proteinsA$ecoliNoTM[3]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff200>0 & ecoliA$TM==0))
proteinsA$ecoli[4]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff100>0 & ecoliA$TM==0))

proteinsA$sacch=NULL
proteinsA$sacchTM[1]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff1000>0 & sacchA$TM==1))
proteinsA$sacchTM[2]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff500>0 & sacchA$TM==1))
proteinsA$sacchTM[3]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff200>0 & sacchA$TM==1))
proteinsA$sacchTM[4]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff100>0 & sacchA$TM==1))
proteinsA$sacchNoTM[1]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff1000>0 & sacchA$TM==0))
proteinsA$sacchNoTM[2]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff500>0 & sacchA$TM==0))
proteinsA$sacchNoTM[3]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff200>0 & sacchA$TM==0))
proteinsA$sacch[4]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff100>0 & sacchA$TM==0))

proteinsA$homo=NULL
proteinsA$homoTM[1]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff1000>0 & homoA$TM==1))
proteinsA$homoTM[2]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff500>0 & homoA$TM==1))
proteinsA$homoTM[3]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff200>0 & homoA$TM==1))
proteinsA$homoTM[4]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff100>0 & homoA$TM==1))
proteinsA$homoNoTM[1]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff1000>0 & homoA$TM==0))
proteinsA$homoNoTM[2]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff500>0 & homoA$TM==0))
proteinsA$homoNoTM[3]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff200>0 & homoA$TM==0))
proteinsA$homo[4]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff100>0 & homoA$TM==0))


prA=NULL
prA$Meff1000[1]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff1000>0))
prA$Meff500[1]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff500>0 ))
prA$Meff200[1]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff200>0 ))
prA$Meff100[1]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff100>0 ))
prA$Meff1000[2]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff1000>0))
prA$Meff500[2]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff500>0 ))
prA$Meff200[2]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff200>0 ))
prA$Meff100[2]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff100>0 ))
prA$Meff1000[3]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff1000>0))
prA$Meff500[3]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff500>0 ))
prA$Meff200[3]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff200>0 ))
prA$Meff100[3]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff100>0 ))

# This is proteins without PDB and for different Meff
proteinsA$ecoli[1]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff1000>0))
proteinsA$ecoli[2]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff500>0 ))
proteinsA$ecoli[3]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff200>0 ))
proteinsA$ecoli[4]=length(which(ecoliA$PDB_count == 0 & ecoliA$Pfam_Meff100>0 ))
proteinsA$sacch[1]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff1000>0))
proteinsA$sacch[2]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff500>0 ))
proteinsA$sacch[3]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff200>0 ))
proteinsA$sacch[4]=length(which(sacchA$PDB_count == 0 & sacchA$Pfam_Meff100>0 ))
proteinsA$homo[1]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff1000>0))
proteinsA$homo[2]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff500>0 ))
proteinsA$homo[3]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff200>0 ))
proteinsA$homo[4]=length(which(homoA$PDB_count == 0 & homoA$Pfam_Meff100>0 ))



proteinsB=NULL
proteinsB$count=NULL
proteinsB$count[1]=1000
proteinsB$count[2]=500
proteinsB$count[3]=200
proteinsB$count[4]=100
proteinsB$ecoli=NULL
proteinsB$ecoliTM[1]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff1000>0 & ecoliB$TM==1))
proteinsB$ecoliTM[2]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff500>0 & ecoliB$TM==1))
proteinsB$ecoliTM[3]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff200>0 & ecoliB$TM==1))
proteinsB$ecoliTM[4]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff100>0 & ecoliB$TM==1))
proteinsB$ecoliNoTM[1]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff1000>0 & ecoliB$TM==0))
proteinsB$ecoliNoTM[2]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff500>0 & ecoliB$TM==0))
proteinsB$ecoliNoTM[3]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff200>0 & ecoliB$TM==0))
proteinsB$ecoli[4]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff100>0 & ecoliB$TM==0))

proteinsB$sacch=NULL
proteinsB$sacchTM[1]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff1000>0 & sacchB$TM==1))
proteinsB$sacchTM[2]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff500>0 & sacchB$TM==1))
proteinsB$sacchTM[3]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff200>0 & sacchB$TM==1))
proteinsB$sacchTM[4]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff100>0 & sacchB$TM==1))
proteinsB$sacchNoTM[1]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff1000>0 & sacchB$TM==0))
proteinsB$sacchNoTM[2]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff500>0 & sacchB$TM==0))
proteinsB$sacchNoTM[3]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff200>0 & sacchB$TM==0))
proteinsB$sacch[4]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff100>0 & sacchB$TM==0))

proteinsB$homo=NULL
proteinsB$homoTM[1]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff1000>0 & homoB$TM==1))
proteinsB$homoTM[2]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff500>0 & homoB$TM==1))
proteinsB$homoTM[3]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff200>0 & homoB$TM==1))
proteinsB$homoTM[4]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff100>0 & homoB$TM==1))
proteinsB$homoNoTM[1]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff1000>0 & homoB$TM==0))
proteinsB$homoNoTM[2]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff500>0 & homoB$TM==0))
proteinsB$homoNoTM[3]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff200>0 & homoB$TM==0))
proteinsB$homo[4]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff100>0 & homoB$TM==0))


prB=NULL
prB$Meff1000[1]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff1000>0))
prB$Meff500[1]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff500>0 ))
prB$Meff200[1]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff200>0 ))
prB$Meff100[1]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff100>0 ))
prB$Meff1000[2]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff1000>0))
prB$Meff500[2]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff500>0 ))
prB$Meff200[2]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff200>0 ))
prB$Meff100[2]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff100>0 ))
prB$Meff1000[3]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff1000>0))
prB$Meff500[3]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff500>0 ))
prB$Meff200[3]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff200>0 ))
prB$Meff100[3]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff100>0 ))

proteinsB$ecoli[1]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff1000>0))
proteinsB$ecoli[2]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff500>0 ))
proteinsB$ecoli[3]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff200>0 ))
proteinsB$ecoli[4]=length(which(ecoliB$PDB_count == 0 & ecoliB$Pfam_Meff100>0 ))
proteinsB$sacch[1]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff1000>0))
proteinsB$sacch[2]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff500>0 ))
proteinsB$sacch[3]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff200>0 ))
proteinsB$sacch[4]=length(which(sacchB$PDB_count == 0 & sacchB$Pfam_Meff100>0 ))
proteinsB$homo[1]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff1000>0))
proteinsB$homo[2]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff500>0 ))
proteinsB$homo[3]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff200>0 ))
proteinsB$homo[4]=length(which(homoB$PDB_count == 0 & homoB$Pfam_Meff100>0 ))



proteinsC=NULL
proteinsC$count=NULL
proteinsC$count[1]=1000
proteinsC$count[2]=500
proteinsC$count[3]=200
proteinsC$count[4]=100
proteinsC$ecoli=NULL
proteinsC$ecoliTM[1]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff1000>0 & ecoliC$TM==1))
proteinsC$ecoliTM[2]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff500>0 & ecoliC$TM==1))
proteinsC$ecoliTM[3]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff200>0 & ecoliC$TM==1))
proteinsC$ecoliTM[4]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff100>0 & ecoliC$TM==1))
proteinsC$ecoliNoTM[1]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff1000>0 & ecoliC$TM==0))
proteinsC$ecoliNoTM[2]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff500>0 & ecoliC$TM==0))
proteinsC$ecoliNoTM[3]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff200>0 & ecoliC$TM==0))
proteinsC$ecoli[4]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff100>0 & ecoliC$TM==0))

proteinsC$sacch=NULL
proteinsC$sacchTM[1]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff1000>0 & sacchC$TM==1))
proteinsC$sacchTM[2]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff500>0 & sacchC$TM==1))
proteinsC$sacchTM[3]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff200>0 & sacchC$TM==1))
proteinsC$sacchTM[4]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff100>0 & sacchC$TM==1))
proteinsC$sacchNoTM[1]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff1000>0 & sacchC$TM==0))
proteinsC$sacchNoTM[2]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff500>0 & sacchC$TM==0))
proteinsC$sacchNoTM[3]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff200>0 & sacchC$TM==0))
proteinsC$sacch[4]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff100>0 & sacchC$TM==0))

proteinsC$homo=NULL
proteinsC$homoTM[1]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff1000>0 & homoC$TM==1))
proteinsC$homoTM[2]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff500>0 & homoC$TM==1))
proteinsC$homoTM[3]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff200>0 & homoC$TM==1))
proteinsC$homoTM[4]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff100>0 & homoC$TM==1))
proteinsC$homoNoTM[1]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff1000>0 & homoC$TM==0))
proteinsC$homoNoTM[2]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff500>0 & homoC$TM==0))
proteinsC$homoNoTM[3]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff200>0 & homoC$TM==0))
proteinsC$homo[4]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff100>0 & homoC$TM==0))


prC=NULL
prC$Meff1000[1]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff1000>0))
prC$Meff500[1]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff500>0 ))
prC$Meff200[1]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff200>0 ))
prC$Meff100[1]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff100>0 ))
prC$Meff1000[2]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff1000>0))
prC$Meff500[2]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff500>0 ))
prC$Meff200[2]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff200>0 ))
prC$Meff100[2]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff100>0 ))
prC$Meff1000[3]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff1000>0))
prC$Meff500[3]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff500>0 ))
prC$Meff200[3]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff200>0 ))
prC$Meff100[3]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff100>0 ))

proteinsC$ecoli[1]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff1000>0))
proteinsC$ecoli[2]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff500>0 ))
proteinsC$ecoli[3]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff200>0 ))
proteinsC$ecoli[4]=length(which(ecoliC$PDB_count == 0 & ecoliC$Pfam_Meff100>0 ))
proteinsC$sacch[1]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff1000>0))
proteinsC$sacch[2]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff500>0 ))
proteinsC$sacch[3]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff200>0 ))
proteinsC$sacch[4]=length(which(sacchC$PDB_count == 0 & sacchC$Pfam_Meff100>0 ))
proteinsC$homo[1]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff1000>0))
proteinsC$homo[2]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff500>0 ))
proteinsC$homo[3]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff200>0 ))
proteinsC$homo[4]=length(which(homoC$PDB_count == 0 & homoC$Pfam_Meff100>0 ))




                                        # Calculate fractions
# PDB
EcoliANoModel=sum(ecoliATMAll[1:2],ecoliATMAll[7:8],ecoliATMAll[13:14])/sum(ecoliATMAll)*100
EcoliAPDB=sum(ecoliATMAll[19:24])/sum(ecoliATMAll)*100
EcoliAGLOB=sum(ecoliATMAll[3:6])/sum(ecoliATMAll)*100
EcoliATM=sum(ecoliATMAll[9:12])/sum(ecoliATMAll)*100
EcoliATMall=sum(ecoliATMAll[7:12])/sum(ecoliATMAll)*100
EcoliADISO=sum(ecoliATMAll[15:18])/sum(ecoliATMAll)*100
EcoliADISOall=sum(ecoliATMAll[13:18])/sum(ecoliATMAll)*100
EcoliAHundred=sum(ecoliATMAll[3:6],ecoliATMAll[9:12],ecoliATMAll[15:18])/sum(ecoliATMAll)*100
EcoliAThousand=sum(ecoliATMAll[6],ecoliATMAll[12],ecoliATMAll[18])/sum(ecoliATMAll)*100
EcoliAPDB
EcoliAHundred
EcoliAPDB+EcoliAHundred

                                        # PDB
SacchANoModel=sum(sacchATMAll[1:2],sacchATMAll[7:8],sacchATMAll[13:14])/sum(sacchATMAll)*100
SacchAPDB=sum(sacchATMAll[19:24])/sum(sacchATMAll)*100
SacchAGLOB=sum(sacchATMAll[3:6])/sum(sacchATMAll)*100
SacchATM=sum(sacchATMAll[9:12])/sum(sacchATMAll)*100
SacchATMall=sum(sacchATMAll[7:12])/sum(sacchATMAll)*100
SacchADISO=sum(sacchATMAll[15:18])/sum(sacchATMAll)*100
SacchADISOall=sum(sacchATMAll[13:18])/sum(sacchATMAll)*100
SacchAHundred=sum(sacchATMAll[3:6],sacchATMAll[9:12],sacchATMAll[15:18])/sum(sacchATMAll)*100
SacchAThousand=sum(sacchATMAll[6],sacchATMAll[12],sacchATMAll[18])/sum(sacchATMAll)*100
SacchAPDB
SacchAHundred
SacchAPDB+SacchAHundred
                                        # PDB
HomoANoModel=sum(homoATMAll[1:2],homoATMAll[7:8],homoATMAll[13:14])/sum(homoATMAll)*100
HomoAPDB=sum(homoATMAll[19:24])/sum(homoATMAll)*100
HomoAGLOB=sum(homoATMAll[3:6])/sum(homoATMAll)*100
HomoATM=sum(homoATMAll[9:12])/sum(homoATMAll)*100
HomoATMall=sum(homoATMAll[7:12])/sum(homoATMAll)*100
HomoADISO=sum(homoATMAll[15:18])/sum(homoATMAll)*100
HomoADISOall=sum(homoATMAll[13:18])/sum(homoATMAll)*100
HomoAHundred=sum(homoATMAll[3:6],homoATMAll[9:12],homoATMAll[15:18])/sum(homoATMAll)*100
HomoAThousand=sum(homoATMAll[6],homoATMAll[12],homoATMAll[18])/sum(homoATMAll)*100
HomoAPDB
HomoAHundred
HomoAPDB+HomoAHundred



                                        # Calculate fractions
# PDB
EcoliBNoModel=sum(ecoliBTMAll[1:2],ecoliBTMAll[7:8],ecoliBTMAll[13:14])/sum(ecoliBTMAll)*100
EcoliBPDB=sum(ecoliBTMAll[19:24])/sum(ecoliBTMAll)*100
EcoliBGLOB=sum(ecoliBTMAll[3:6])/sum(ecoliBTMAll)*100
EcoliBTM=sum(ecoliBTMAll[9:12])/sum(ecoliBTMAll)*100
EcoliBTMall=sum(ecoliBTMAll[7:12])/sum(ecoliBTMAll)*100
EcoliBDISO=sum(ecoliBTMAll[15:18])/sum(ecoliBTMAll)*100
EcoliBDISOall=sum(ecoliBTMAll[13:18])/sum(ecoliBTMAll)*100
EcoliBHundred=sum(ecoliBTMAll[3:6],ecoliBTMAll[9:12],ecoliBTMAll[15:18])/sum(ecoliBTMAll)*100
EcoliBThousand=sum(ecoliBTMAll[6],ecoliBTMAll[12],ecoliBTMAll[18])/sum(ecoliBTMAll)*100
EcoliBPDB
EcoliBHundred
EcoliBPDB+EcoliBHundred

                                        # PDB
SacchBNoModel=sum(sacchBTMAll[1:2],sacchBTMAll[7:8],sacchBTMAll[13:14])/sum(sacchBTMAll)*100
SacchBPDB=sum(sacchBTMAll[19:24])/sum(sacchBTMAll)*100
SacchBGLOB=sum(sacchBTMAll[3:6])/sum(sacchBTMAll)*100
SacchBTM=sum(sacchBTMAll[9:12])/sum(sacchBTMAll)*100
SacchBTMall=sum(sacchBTMAll[7:12])/sum(sacchBTMAll)*100
SacchBDISO=sum(sacchBTMAll[15:18])/sum(sacchBTMAll)*100
SacchBDISOall=sum(sacchBTMAll[13:18])/sum(sacchBTMAll)*100
SacchBHundred=sum(sacchBTMAll[3:6],sacchBTMAll[9:12],sacchBTMAll[15:18])/sum(sacchBTMAll)*100
SacchBThousand=sum(sacchBTMAll[6],sacchBTMAll[12],sacchBTMAll[18])/sum(sacchBTMAll)*100
SacchBPDB
SacchBHundred
SacchBPDB+SacchBHundred
                                        # PDB
HomoBNoModel=sum(homoBTMAll[1:2],homoBTMAll[7:8],homoBTMAll[13:14])/sum(homoBTMAll)*100
HomoBPDB=sum(homoBTMAll[19:24])/sum(homoBTMAll)*100
HomoBGLOB=sum(homoBTMAll[3:6])/sum(homoBTMAll)*100
HomoBTM=sum(homoBTMAll[9:12])/sum(homoBTMAll)*100
HomoBTMall=sum(homoBTMAll[7:12])/sum(homoBTMAll)*100
HomoBDISO=sum(homoBTMAll[15:18])/sum(homoBTMAll)*100
HomoBDISOall=sum(homoBTMAll[13:18])/sum(homoBTMAll)*100
HomoBHundred=sum(homoBTMAll[3:6],homoBTMAll[9:12],homoBTMAll[15:18])/sum(homoBTMAll)*100
HomoBThousand=sum(homoBTMAll[6],homoBTMAll[12],homoBTMAll[18])/sum(homoBTMAll)*100
HomoBPDB
HomoBHundred
HomoBPDB+HomoBHundred


                                        # Calculate fractions
# PDB
EcoliCNoModel=sum(ecoliCTMAll[1:2],ecoliCTMAll[7:8],ecoliCTMAll[13:14])/sum(ecoliCTMAll)*100
EcoliCPDB=sum(ecoliCTMAll[19:24])/sum(ecoliCTMAll)*100
EcoliCGLOB=sum(ecoliCTMAll[3:6])/sum(ecoliCTMAll)*100
EcoliCTM=sum(ecoliCTMAll[9:12])/sum(ecoliCTMAll)*100
EcoliCTMall=sum(ecoliCTMAll[7:12])/sum(ecoliCTMAll)*100
EcoliCDISO=sum(ecoliCTMAll[15:18])/sum(ecoliCTMAll)*100
EcoliCDISOall=sum(ecoliCTMAll[13:18])/sum(ecoliCTMAll)*100
EcoliCHundred=sum(ecoliCTMAll[3:6],ecoliCTMAll[9:12],ecoliCTMAll[15:18])/sum(ecoliCTMAll)*100
EcoliCThousand=sum(ecoliCTMAll[6],ecoliCTMAll[12],ecoliCTMAll[18])/sum(ecoliCTMAll)*100
EcoliCPDB
EcoliCHundred
EcoliCPDB+EcoliCHundred

                                        # PDB
SacchCNoModel=sum(sacchCTMAll[1:2],sacchCTMAll[7:8],sacchCTMAll[13:14])/sum(sacchCTMAll)*100
SacchCPDB=sum(sacchCTMAll[19:24])/sum(sacchCTMAll)*100
SacchCGLOB=sum(sacchCTMAll[3:6])/sum(sacchCTMAll)*100
SacchCTM=sum(sacchCTMAll[9:12])/sum(sacchCTMAll)*100
SacchCTMall=sum(sacchCTMAll[7:12])/sum(sacchCTMAll)*100
SacchCDISO=sum(sacchCTMAll[15:18])/sum(sacchCTMAll)*100
SacchCDISOall=sum(sacchCTMAll[13:18])/sum(sacchCTMAll)*100
SacchCHundred=sum(sacchCTMAll[3:6],sacchCTMAll[9:12],sacchCTMAll[15:18])/sum(sacchCTMAll)*100
SacchCThousand=sum(sacchCTMAll[6],sacchCTMAll[12],sacchCTMAll[18])/sum(sacchCTMAll)*100
SacchCPDB
SacchCHundred
SacchCPDB+SacchCHundred
                                        # PDB
HomoCNoModel=sum(homoCTMAll[1:2],homoCTMAll[7:8],homoCTMAll[13:14])/sum(homoCTMAll)*100
HomoCPDB=sum(homoCTMAll[19:24])/sum(homoCTMAll)*100
HomoCGLOB=sum(homoCTMAll[3:6])/sum(homoCTMAll)*100
HomoCTM=sum(homoCTMAll[9:12])/sum(homoCTMAll)*100
HomoCTMall=sum(homoCTMAll[7:12])/sum(homoCTMAll)*100
HomoCDISO=sum(homoCTMAll[15:18])/sum(homoCTMAll)*100
HomoCDISOall=sum(homoCTMAll[13:18])/sum(homoCTMAll)*100
HomoCHundred=sum(homoCTMAll[3:6],homoCTMAll[9:12],homoCTMAll[15:18])/sum(homoCTMAll)*100
HomoCThousand=sum(homoCTMAll[6],homoCTMAll[12],homoCTMAll[18])/sum(homoCTMAll)*100
HomoCPDB
HomoCHundred
HomoCPDB+HomoCHundred


labelssort=NULL
labelssort<-append(labelssort,labels[1:2])
labelssort<-append(labelssort,labels[7:8])
labelssort<-append(labelssort,labels[13:14])
labelssort<-append(labelssort,labels[3:6])
labelssort<-append(labelssort,labels[9:12])
labelssort<-append(labelssort,labels[15:18])
labelssort<-append(labelssort,labels[19:24])

colorssort=NULL
colorssort<-append(colorssort,colors[1:2])
colorssort<-append(colorssort,colors[7:8])
colorssort<-append(colorssort,colors[13:14])
colorssort<-append(colorssort,colors[3:6])
colorssort<-append(colorssort,colors[9:12])
colorssort<-append(colorssort,colors[15:18])
colorssort<-append(colorssort,colors[19:24])


ecoliAsort=NULL
ecoliAsort<-append(ecoliAsort,ecoliATMAll[1:2])
ecoliAsort<-append(ecoliAsort,ecoliATMAll[7:8])
ecoliAsort<-append(ecoliAsort,ecoliATMAll[13:14])
ecoliAsort<-append(ecoliAsort,ecoliATMAll[3:6])
ecoliAsort<-append(ecoliAsort,ecoliATMAll[9:12])
ecoliAsort<-append(ecoliAsort,ecoliATMAll[15:18])
ecoliAsort<-append(ecoliAsort,ecoliATMAll[19:24])

ecoliBsort=NULL
ecoliBsort<-append(ecoliBsort,ecoliBTMAll[1:2])
ecoliBsort<-append(ecoliBsort,ecoliBTMAll[7:8])
ecoliBsort<-append(ecoliBsort,ecoliBTMAll[13:14])
ecoliBsort<-append(ecoliBsort,ecoliBTMAll[3:6])
ecoliBsort<-append(ecoliBsort,ecoliBTMAll[9:12])
ecoliBsort<-append(ecoliBsort,ecoliBTMAll[15:18])
ecoliBsort<-append(ecoliBsort,ecoliBTMAll[19:24])


ecoliCsort=NULL
ecoliCsort<-append(ecoliCsort,ecoliCTMAll[1:2])
ecoliCsort<-append(ecoliCsort,ecoliCTMAll[7:8])
ecoliCsort<-append(ecoliCsort,ecoliCTMAll[13:14])
ecoliCsort<-append(ecoliCsort,ecoliCTMAll[3:6])
ecoliCsort<-append(ecoliCsort,ecoliCTMAll[9:12])
ecoliCsort<-append(ecoliCsort,ecoliCTMAll[15:18])
ecoliCsort<-append(ecoliCsort,ecoliCTMAll[19:24])



sacchAsort=NULL
sacchAsort<-append(sacchAsort,sacchATMAll[1:2])
sacchAsort<-append(sacchAsort,sacchATMAll[7:8])
sacchAsort<-append(sacchAsort,sacchATMAll[13:14])
sacchAsort<-append(sacchAsort,sacchATMAll[3:6])
sacchAsort<-append(sacchAsort,sacchATMAll[9:12])
sacchAsort<-append(sacchAsort,sacchATMAll[15:18])
sacchAsort<-append(sacchAsort,sacchATMAll[19:24])

sacchBsort=NULL
sacchBsort<-append(sacchBsort,sacchBTMAll[1:2])
sacchBsort<-append(sacchBsort,sacchBTMAll[7:8])
sacchBsort<-append(sacchBsort,sacchBTMAll[13:14])
sacchBsort<-append(sacchBsort,sacchBTMAll[3:6])
sacchBsort<-append(sacchBsort,sacchBTMAll[9:12])
sacchBsort<-append(sacchBsort,sacchBTMAll[15:18])
sacchBsort<-append(sacchBsort,sacchBTMAll[19:24])


sacchCsort=NULL
sacchCsort<-append(sacchCsort,sacchCTMAll[1:2])
sacchCsort<-append(sacchCsort,sacchCTMAll[7:8])
sacchCsort<-append(sacchCsort,sacchCTMAll[13:14])
sacchCsort<-append(sacchCsort,sacchCTMAll[3:6])
sacchCsort<-append(sacchCsort,sacchCTMAll[9:12])
sacchCsort<-append(sacchCsort,sacchCTMAll[15:18])
sacchCsort<-append(sacchCsort,sacchCTMAll[19:24])



homoAsort=NULL
homoAsort<-append(homoAsort,homoATMAll[1:2])
homoAsort<-append(homoAsort,homoATMAll[7:8])
homoAsort<-append(homoAsort,homoATMAll[13:14])
homoAsort<-append(homoAsort,homoATMAll[3:6])
homoAsort<-append(homoAsort,homoATMAll[9:12])
homoAsort<-append(homoAsort,homoATMAll[15:18])
homoAsort<-append(homoAsort,homoATMAll[19:24])

homoBsort=NULL
homoBsort<-append(homoBsort,homoBTMAll[1:2])
homoBsort<-append(homoBsort,homoBTMAll[7:8])
homoBsort<-append(homoBsort,homoBTMAll[13:14])
homoBsort<-append(homoBsort,homoBTMAll[3:6])
homoBsort<-append(homoBsort,homoBTMAll[9:12])
homoBsort<-append(homoBsort,homoBTMAll[15:18])
homoBsort<-append(homoBsort,homoBTMAll[19:24])


homoCsort=NULL
homoCsort<-append(homoCsort,homoCTMAll[1:2])
homoCsort<-append(homoCsort,homoCTMAll[7:8])
homoCsort<-append(homoCsort,homoCTMAll[13:14])
homoCsort<-append(homoCsort,homoCTMAll[3:6])
homoCsort<-append(homoCsort,homoCTMAll[9:12])
homoCsort<-append(homoCsort,homoCTMAll[15:18])
homoCsort<-append(homoCsort,homoCTMAll[19:24])



genome="Ecoli-prot-by-year"
outfile=paste("figures/",genome,"-Detailed-pie.png",sep="")
pct <- round(ecoliCTMAll/sum(ecoliCTMAll)*100,digits=1)
pct <- round(sacchCDisoAll/sum(sacchCDisoAll)*100,digits=1)
Names <- paste(labelssort,pct)
names <- paste(Names,"%",sep="")
png(outfile,width=1280,height=1280)
pie(ecoliCsort, labels=names,col=colorssort,sub="Fraction of proteins",main="E Coli",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0,cex.lab=0.4)
floating.pie(0,0,ecoliBsort, col=colorssort,main='',radius=0.7)
floating.pie(0,0,ecoliAsort, col=colorssort,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"2015",cex=2.0,border=NA)
legend(-.7,.0,"2010",cex=2.0,border=NA)
legend(-.4,.0,"2005",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Detailed-bar.png",sep="")
png(outfile,width=1280,height=1280)
test=matrix(c(rev(ecoliAsort),rev(ecoliBsort),rev(ecoliCsort)),nrow=24,ncol=3)
barplot(test,col=rev(colorssort),main="E. Coli",legend=rev(names),xlim=c(0,4.5),xlab="Year",ylab="Number or proteins",names=years,cex.names=2,cex.axis=2.,cex=2)
dev.off()

genome="Yeast-prot-by-year"
outfile=paste("figures/",genome,"-Detailed-pie.png",sep="")
pct <- round(sacchCDisoAll/sum(sacchCDisoAll)*100,digits=1)
Names <- paste(labelssort,pct)
names <- paste(Names,"%",sep="")
png(outfile,width=1280,height=1280)
pie(sacchCsort, labels=labelssort,col=colorssort,sub="Fraction of proteins",main="Yeast",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,sacchBsort, col=colorssort,main='',radius=0.7)
floating.pie(0,0,sacchAsort, col=colorssort,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"2015",cex=2.0,border=NA)
legend(-.7,.0,"2010",cex=2.0,border=NA)
legend(-.4,.0,"2005",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Detailed-bar.png",sep="")
png(outfile,width=1280,height=1280)
test=matrix(c(rev(sacchAsort),rev(sacchBsort),rev(sacchCsort)),nrow=24,ncol=3)
barplot(test,col=rev(colorssort),main="Yeast",legend=rev(names),xlim=c(0,4.5),xlab="Year",ylab="Number or proteins",names=years,cex.names=2,cex.axis=2.,cex=2)
dev.off()

genome="Homo-prot-by-year"
outfile=paste("figures/",genome,"-Detailed-pie.png",sep="")
pct <- round(homoCsort/sum(homoCsort)*100,digits=1)
Names <- paste(labelssort,pct)
names <- paste(Names,"%",sep="")
png(outfile,width=1280,height=1280)
pie(homoCsort, labels=names,col=colorssort,sub="Fraction of proteins",main="Homo Sapiens",radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
floating.pie(0,0,homoBsort, col=colorssort,main='',radius=0.7)
floating.pie(0,0,homoAsort, col=colorssort,main='',radius=0.4)
floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
legend(-1,.0,"2015",cex=2.0,border=NA)
legend(-.7,.0,"2010",cex=2.0,border=NA)
legend(-.4,.0,"2005",cex=2.0,border=NA)
dev.off()

outfile=paste("figures/",genome,"-Detailed-bar.png",sep="")
png(outfile,width=1280,height=1280)
test=matrix(c(rev(homoAsort),rev(homoBsort),rev(homoCsort)),nrow=24,ncol=3)
barplot(test,col=rev(colorssort),main="Homo Sapiens",legend=rev(names),xlim=c(0,4.5),xlab="Year",ylab="Number or proteins",names=years,cex.names=2,cex.axis=2.,cex=2)
dev.off()


genome="All"
#outfile=paste("figures/",genome,"-2015-Detailed-bar.png",sep="")
outfile=paste("figures/",genome,"-2015-Detailed-bar.eps",sep="")
png(outfile,width=1280,height=1280)
genlabels=NULL
genlabels[1]="Homo Sapiens"
genlabels[2]="Yeast"
genlabels[3]="E.coli"
ecoliCsortNorm <- ecoliCsort/sum(ecoliCsort)
sacchCsortNorm <- sacchCsort/sum(sacchCsort)
homoCsortNorm <- homoCsort/sum(homoCsort)
test=matrix(c(rev(homoCsortNorm),rev(sacchCsortNorm),rev(ecoliCsortNorm)),nrow=24,ncol=3)
bp<-barplot(test,col=rev(colorssort),main="Genomes 2015",legend=rev(labelssort),xlim=c(0,4.5),xlab="Year",ylab="Number or proteins",names=genlabels,cex.names=2,cex.axis=2.,cex=2)

 

fractionCPDB=NULL
fractionCPDB[1]=round(HomoCPDB,0)
fractionCPDB[2]=round(SacchCPDB,0)
fractionCPDB[3]=round(EcoliCPDB,0)
labelCPDB <- paste(fractionCPDB,"%",sep="")
text(bp, 0, labelCPDB,cex=2,pos=3,col="white") 



fractionCDISO=NULL
fractionCDISOall=NULL
fractionCDISO[1]=round(HomoCDISO,1)
fractionCDISO[2]=round(SacchCDISO,1)
fractionCDISO[3]=round(EcoliCDISO,1)
fractionCDISOall[1]=round(HomoCDISOall,1)
fractionCDISOall[2]=round(SacchCDISOall,1)
fractionCDISOall[3]=round(EcoliCDISOall,1)
labelCDISO <- paste(fractionCDISO,"%",sep="")
text(bp, fractionCPDB/100 , labelCDISO,cex=2,pos=3,col="black") 


fractionCTM=NULL
fractionCTMall=NULL
fractionCTM[1]=round(HomoCTM,0)
fractionCTM[2]=round(SacchCTM,0)
fractionCTM[3]=round(EcoliCTM,0)
fractionCTMall[1]=round(HomoCTMall,1)
fractionCTMall[2]=round(SacchCTMall,1)
fractionCTMall[3]=round(EcoliCTMall,1)
labelCTM <- paste(fractionCTM,"%",sep="")
text(bp, (fractionCPDB+fractionCDISO )/100, labelCTM,cex=2,pos=3,col="yellow") 

fractionCGLOB=NULL
fractionCGLOBall=NULL
fractionCGLOB[1]=round(HomoCGLOB,0)
fractionCGLOB[2]=round(SacchCGLOB,0)
fractionCGLOB[3]=round(EcoliCGLOB,0)
labelCGLOB <- paste(fractionCGLOB,"%",sep="")
text(bp, ( fractionCPDB+fractionCDISO+fractionCTM )/100,labelCGLOB,cex=2,pos=3,col="red") 

dev.off()
                                        #---------New start
