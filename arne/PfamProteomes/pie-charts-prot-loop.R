library(vioplot)
library('plotrix')
genomes=NULL
years=NULL
years[1]="2005"
years[2]="2010"
years[3]="2015"
genomes[1]="escherichia_coli"
genomes[2]="saccharomyces_cerevisae"
genomes[3]="homo_sapiens"

for (year in years){

    file=paste("data/",genomes[1],".df.",year,".blast.prot.tsv",sep="")
    ecoli<-read.table(file, sep='\t', header=T)
    file=paste("data/",genomes[2],".df.",year,".blast.prot.tsv",sep="")
    sacch<-read.table(file, sep='\t', header=T)
    file=paste("data/",genomes[3],".df.",year,".blast.prot.tsv",sep="")
    homo<-read.table(file, sep='\t', header=T)


    columns=NULL
    cutoffs=NULL
    
    cutoffs[1]=-1
    cutoffs[2]=0
    cutoffs[3]=100
    cutoffs[4]=200
    cutoffs[5]=500
    cutoffs[6]=1000
    cutoffs[7]=1000000000
    
# ----- Pre processing ---

    ecoli$Pfam_Meff=ecoli$Pfam_Meff100
    ecoli$Pfam_Meff[ecoli$Pfam_Meff100>0]<-101
    ecoli$Pfam_Meff[ecoli$Pfam_Meff200>0]<-201
    ecoli$Pfam_Meff[ecoli$Pfam_Meff500>0]<-501
    ecoli$Pfam_Meff[ecoli$Pfam_Meff1000>0]<-1001
    ecoli$PDB_ID=ecoli$PDB_count
    ecoli$PDB_ID[ecoli$PDB_count==0]<-""
    ecoli$PDB_ID[ecoli$PDB_count>0]<-"PDB"
    
    sacch$Pfam_Meff=sacch$Pfam_Meff100
    sacch$Pfam_Meff[sacch$Pfam_Meff100>0]<-101
    sacch$Pfam_Meff[sacch$Pfam_Meff200>0]<-201
    sacch$Pfam_Meff[sacch$Pfam_Meff500>0]<-501
    sacch$Pfam_Meff[sacch$Pfam_Meff1000>0]<-1001
    sacch$PDB_ID=sacch$PDB_count
    sacch$PDB_ID[sacch$PDB_count==0]<-""
    sacch$PDB_ID[sacch$PDB_count>0]<-"PDB"
    
    homo$Pfam_Meff=homo$Pfam_Meff100
    homo$Pfam_Meff[homo$Pfam_Meff100>0]<-101
    homo$Pfam_Meff[homo$Pfam_Meff200>0]<-201
    homo$Pfam_Meff[homo$Pfam_Meff500>0]<-501
    homo$Pfam_Meff[homo$Pfam_Meff1000>0]<-1001
    homo$PDB_ID=homo$PDB_count
    homo$PDB_ID[homo$PDB_count==0]<-""
    homo$PDB_ID[homo$PDB_count>0]<-"PDB"
    
    
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
    genome="ecoli-prot"
    outfile=paste(genome,".",year,"-PDB-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    fraction=(sum(ecoliPDB)/(sum(ecoliNoPDB)+sum(ecoliPDB)))**2
    pie(ecoliAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
    floating.pie(0,0,ecoliPDB, col=colors,main='',radius=0.7)
    floating.pie(0,0,ecoliNoPDB, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"PDB",cex=2.0,border=NA)
    legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-TM-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(ecoliAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
    floating.pie(0,0,ecoliTM, col=colors,main='',radius=0.7)
    floating.pie(0,0,ecoliNoTM, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"TM",cex=2.0,border=NA)
    legend(-.4,.0,"NoTM",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-Diso1-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(ecoliAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
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
    
    genome="ecoli-all-prot"
    pct <- round(ecoliPDBAll/sum(ecoliAll)*100,digits=1)
    names <- paste(labels,pct)
    names <- paste(names,"%",sep="")
    
    outfile=paste(genome,".",year,"-TM-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(ecoliPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
    floating.pie(0,0,ecoliPDBTM, col=colors,main='',radius=0.7)
    floating.pie(0,0,ecoliPDBNoTM, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"TM",cex=2.0,border=NA)
    legend(-.4,.0,"NoTM",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-Diso2-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(ecoliPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
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
    genome="sacch-prot"
    outfile=paste(genome,".",year,"-PDB.png",sep="")
    png(outfile,width=1280,height=1280)
    fraction=(sum(sacchPDB)/(sum(sacchNoPDB)+sum(sacchPDB)))**2
    pie(sacchAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
    floating.pie(0,0,sacchPDB, col=colors,main='',radius=0.7)
    floating.pie(0,0,sacchNoPDB, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"PDB",cex=2.0,border=NA)
    legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-TM-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(sacchAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
    floating.pie(0,0,sacchTM, col=colors,main='',radius=0.7)
    floating.pie(0,0,sacchNoTM, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"TM",cex=2.0,border=NA)
    legend(-.4,.0,"NoTM",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-Diso1-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(sacchAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
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
    
    genome="sacch-all-prot"
    pct <- round(sacchAll/sum(sacchPDBAll)*100,digits=1)
    names <- paste(labels,pct)
    names <- paste(names,"%",sep="")
    
    outfile=paste(genome,".",year,"-TM-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(sacchPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
    floating.pie(0,0,sacchPDBTM, col=colors,main='',radius=0.7)
    floating.pie(0,0,sacchPDBNoTM, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"TM",cex=2.0,border=NA)
    legend(-.4,.0,"NoTM",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-Diso2-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(sacchPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
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
    genome="homo-prot"
    outfile=paste(genome,".",year,"-PDB.png",sep="")
    png(outfile,width=1280,height=1280)
    fraction=(sum(homoPDB)/(sum(homoNoPDB)+sum(homoPDB)))**2
    pie(homoAll, labels=names,col=colors,main=genome,sub="PDB",radius=iniR,border = NA)
    floating.pie(0,0,homoPDB, col=colors,main='',radius=0.7)
    floating.pie(0,0,homoNoPDB, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"PDB",cex=2.0,border=NA)
    legend(-.4,.0,"NoPDB",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-TM-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(homoAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
    floating.pie(0,0,homoTM, col=colors,main='',radius=0.7)
    floating.pie(0,0,homoNoTM, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"TM",cex=2.0,border=NA)
    legend(-.4,.0,"NoTM",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-Diso1-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(homoAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
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
    
    genome="homo-all-prot"
    pct <- round(homoAll/sum(homoPDBAll)*100,digits=1)
    names <- paste(labels,pct)
    names <- paste(names,"%",sep="")
    
    outfile=paste(genome,".",year,"-TM-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(homoPDBAll, labels=names,col=colors,main=genome,sub="TM",radius=iniR,border = NA)
    floating.pie(0,0,homoPDBTM, col=colors,main='',radius=0.7)
    floating.pie(0,0,homoPDBNoTM, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"TM",cex=2.0,border=NA)
    legend(-.4,.0,"NoTM",cex=2.0,border=NA)
    dev.off()
    
    outfile=paste(genome,".",year,"-Diso2-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    pie(homoPDBAll, labels=names,col=colors,main=genome,sub="Diso",radius=iniR,border = NA)
    floating.pie(0,0,homoPDBDiso, col=colors,main='',radius=0.7)
    floating.pie(0,0,homoPDBNoDiso, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
                                        #legend(0, 0, gsub("_"," ",names(colors)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
    legend(-1,.0,"all",cex=2.0,border=NA)
    legend(-.7,.0,"Diso",cex=2.0,border=NA)
    legend(-.4,.0,"NoDiso",cex=2.0,border=NA)
    dev.off()
    
    
    
                                        #Summary of all
    genome="All-prot"
    outfile=paste(genome,".",year,"-summary-pie.png",sep="")
    png(outfile,width=1280,height=1280)
    
    pie(homoPDBAll, labels=labels,col=colors,sub="Genome Classification",main=year,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
    floating.pie(0,0,sacchPDBAll, col=colors,main='',radius=0.7)
    floating.pie(0,0,ecoliPDBAll, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
    legend(-1,.0,"Homo",cex=2.0,border=NA)
    legend(-.7,.0,"Sacc",cex=2.0,border=NA)
    legend(-.4,.0,"EColi",cex=2.0,border=NA)
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
    
    genome="All-prot"
    outfile=paste(genome,".",year,"-Detailed1-pie.png",sep="")
    
    pct <- round(ecoliDisoAll/sum(ecoliDisoAll)*100,digits=1)
    Names <- paste(labels,pct)
    names <- paste(names,"%",sep="")
    png(outfile,width=1280,height=1280)
    
    pie(homoDisoAll, labels=labels,col=colors,sub="Detailed Genome Classification",main=year,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
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
    
    labels[19]="PDB"
    colors[19]="grey90"
    colors[20]="grey80"
    colors[21]="grey40"
    colors[22]="grey30"
    colors[23]="grey20"
    colors[24]="black"
    
    labels=NULL
    
    labels[1]="Globular"
    labels[7]="TM"
    labels[13]="Disorder"
    labels[19]="PDB"
    
    
    genome="All-prot"
    outfile=paste(genome,".",year,"-Detailed2-pie.png",sep="")
    pct <- round(ecoliDisoAll/sum(ecoliDisoAll)*100,digits=1)
    Names <- paste(labels,pct)
    names <- paste(names,"%",sep="")
    png(outfile,width=1280,height=1280)
    pie(homoTMAll, labels=labels,col=colors,sub="Genome Classification",main=year,radius=iniR,border = NA,cex=2.0,cex.main=3.0,cex.sub=2.0)
    floating.pie(0,0,sacchTMAll, col=colors,main='',radius=0.7)
    floating.pie(0,0,ecoliTMAll, col=colors,main='',radius=0.4)
    floating.pie(0,0,c(1), radius=0.1, col=c('white'), border = NA)
    legend(-1,.0,"Homo",cex=2.0,border=NA)
    legend(-.7,.0,"Sacch",cex=2.0,border=NA)
    legend(-.4,.0,"EColi",cex=2.0,border=NA)
    dev.off()
    
    
    proteins=NULL
    proteins$count=NULL
    proteins$count[1]=1000
    proteins$count[2]=500
    proteins$count[3]=200
    proteins$count[4]=100
    proteins$ecoli=NULL
    proteins$ecoliTM[1]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff1000>0 & ecoli$TM==1))
    proteins$ecoliTM[2]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff500>0 & ecoli$TM==1))
    proteins$ecoliTM[3]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff200>0 & ecoli$TM==1))
    proteins$ecoliTM[4]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff100>0 & ecoli$TM==1))
    proteins$ecoliNoTM[1]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff1000>0 & ecoli$TM==0))
    proteins$ecoliNoTM[2]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff500>0 & ecoli$TM==0))
    proteins$ecoliNoTM[3]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff200>0 & ecoli$TM==0))
    proteins$ecoli[4]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff100>0 & ecoli$TM==0))
    pr=NULL
    pr$Meff1000[1]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff1000>0))
    pr$Meff500[1]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff500>0 ))
    pr$Meff200[1]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff200>0 ))
    pr$Meff100[1]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff100>0 ))
    pr$Meff1000[2]=length(which(sacch$PDB_count == 0 & sacch$Pfam_Meff1000>0))
    pr$Meff500[2]=length(which(sacch$PDB_count == 0 & sacch$Pfam_Meff500>0 ))
    pr$Meff200[2]=length(which(sacch$PDB_count == 0 & sacch$Pfam_Meff200>0 ))
    pr$Meff100[2]=length(which(sacch$PDB_count == 0 & sacch$Pfam_Meff100>0 ))
    pr$Meff1000[3]=length(which(homo$PDB_count == 0 & homo$Pfam_Meff1000>0))
    pr$Meff500[3]=length(which(homo$PDB_count == 0 & homo$Pfam_Meff500>0 ))
    pr$Meff200[3]=length(which(homo$PDB_count == 0 & homo$Pfam_Meff200>0 ))
    pr$Meff100[3]=length(which(homo$PDB_count == 0 & homo$Pfam_Meff100>0 ))
    
    proteins$ecoli[1]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff1000>0))
    proteins$ecoli[2]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff500>0 ))
    proteins$ecoli[3]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff200>0 ))
    proteins$ecoli[4]=length(which(ecoli$PDB_count == 0 & ecoli$Pfam_Meff100>0 ))
    proteins$sacch[1]=length(which(sacch$PDB_count == 0 & sacch$Pfam_Meff1000>0))
    proteins$sacch[2]=length(which(sacch$PDB_count == 0 & sacch$Pfam_Meff500>0 ))
    proteins$sacch[3]=length(which(sacch$PDB_count == 0 & sacch$Pfam_Meff200>0 ))
    proteins$sacch[4]=length(which(sacch$PDB_count == 0 & sacch$Pfam_Meff100>0 ))
    proteins$homo[1]=length(which(homo$PDB_count == 0 & homo$Pfam_Meff1000>0))
    proteins$homo[2]=length(which(homo$PDB_count == 0 & homo$Pfam_Meff500>0 ))
    proteins$homo[3]=length(which(homo$PDB_count == 0 & homo$Pfam_Meff200>0 ))
    proteins$homo[4]=length(which(homo$PDB_count == 0 & homo$Pfam_Meff100>0 ))
    
    genome="Pred-proteins"
    outfile=paste(genome,".",year,"-PDB.png",sep="")
    png(outfile,width=1280,height=1280)
    barplot(c(proteins$ecoli,proteins$sacch,proteins$homo),names.arg=c(proteins$count,proteins$count,proteins$count),sub="Number of proteins that can be predicted for HS, Yeast and E.Coli",main=year,xlab="Meff",ylab="Number or proteins",col=c("darkblue","darkblue","darkblue","darkblue","red","red","red","red","green","green","green","green"),cex.names=2.0,cex.axis=3.0,cex.main=3.0,cex.sub=2.0,ylim=c(0,5000))
    dev.off()
}
    
