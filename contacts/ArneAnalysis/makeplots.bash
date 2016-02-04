#!/bin/bash -x

for i in "$@"
do
    j=`basename $i .PconsC3.l5  `
    j=`basename $j .rr`
    j=`basename $j .gdca`
    j=`basename $j .0.02.plm20`
    f=`echo $i | sed -E "s/.*JH0.001.//"`
    k=`echo $j | sed -E "s/\.fa.*//" | sed -E "s/\.fasta.*//"`
    l=`echo $k | sed -E "s/DisProt-//" | sed -E "s/\-uniprot.*//"`
    if [ -e $k.fa ]
    then
	e="fa"
    else
	e="fasta"
    fi
    if [ -e $l.pdb ]
    then
    ~/git/bioinfo-toolbox/contacts/ArneAnalysis/plot_contact_map.py --cutoff 0.3 --meff $j.gneff --alignment $j.trimmed  --iupred $k-iupred-long.txt  --psipred_vert $j.ss2 --pdb $l.pdb  $k.$e $j.$f
    ~/git/bioinfo-toolbox/contacts/ArneAnalysis/plot_contact_map.py --cutoff 0.3 --alignment $j.trimmed  --iupred $k-iupred-long.txt  --psipred_vert $j.ss2 --pdb $l.pdb  $k.$e $j.$f -o $j-nomeff.pdf
    ~/git/bioinfo-toolbox/contacts/ArneAnalysis/analyse_predictions.py --cutoff 0.3 --meff $j.gneff --alignment $j.trimmed  --iupred $k-iupred-long.txt  --psipred_vert $j.ss2  --pdb $l.pdb $k.$e $j.$f
    else
    ~/git/bioinfo-toolbox/contacts/ArneAnalysis/plot_contact_map.py --cutoff 0.3 --meff $j.gneff --alignment $j.trimmed  --iupred $k-iupred-long.txt  --psipred_vert $j.ss2  $k.$e $j.$f
    ~/git/bioinfo-toolbox/contacts/ArneAnalysis/analyse_predictions.py --cutoff 0.3 --meff $j.gneff --alignment $j.trimmed  --iupred $k-iupred-long.txt  --psipred_vert $j.ss2 $k.$e $j.$f
    fi
done
