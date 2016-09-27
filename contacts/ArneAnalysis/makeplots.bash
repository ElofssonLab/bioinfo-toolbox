#!/bin/bash 

for i in "$@"
do
    j=`basename $i .PconsC3.l5  `
    j=`basename $j .PconsC3.l4  `
    j=`basename $j .PconsC3.l3  `
    j=`basename $j .pconsc3.l5  `
    j=`basename $j .pconsc3.l4  `
    j=`basename $j .pconsc3.l3  `
    j=`basename $j .rr`
    j=`basename $j .gdca`
    j=`basename $j .0.02.plm20`
#    f=`echo $i | sed -E "s/.*JH1.//"`
#    f=`echo $f | sed -E "s/.*HH.//"`
#    f=`echo $i | sed -E "s/.*jhE.//"`
#    f=`echo $f | sed -E "s/.*hhE.//"`
    k=`echo $j | sed -E "s/\.fa.*//" | sed -E "s/\.fasta.*//" | sed -E "s/\.seq.*//"`
    l=`echo $k | sed -E "s/DisProt-//" | sed -E "s/\-uniprot.*//"`
    
    if [ -e $k.fa ]
    then
	e="fa"
    elif [ -e $k.fasta ]
    then
	 e="fasta"
    elif [ -e $k.seq ]
    then
	 e="seq"
    else
	e="fasta"
    fi
    p=' '

    if [ -e $k-iupred-long.txt ]
    then
	p=$p" --iupred $k-iupred-long.txt " 
    fi

#    if [ -e $j.gneff ]
#    then#
#	p=$p" --meff $j.gneff  "
#    fi
    if [ -e $l.pdb ]
    then
	p=$p" --pdb $l.pdb  "
    fi
    if [ -e $j.trimmed ]
    then
	p=$p" --alignment $j.trimmed  "
    fi
    if [ -e $j.ss2 ]
    then
	p=$p" --psipred_vert $j.ss2  "
    fi
    if [ ! -e $i\_ContactMap\.pdf ]
	then
	~/git/bioinfo-toolbox/contacts/ArneAnalysis/plot_contact_map.py --threshold 0.4 --factor 1.0 $p    $k.$e $i
    #~/git/bioinfo-toolbox/contacts/ArneAnalysis/analyse_predictions.py --cutoff 0.4 $p $k.$e $j.$f
    fi
done
