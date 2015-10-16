#!/bin/bash

# Created 2012-02-23, updated 2012-02-23, Nanjiang Shu  
# create simple bar charts dataset for ITOL phylogenetic tree drawing from 
# sorted.orig.topomsa.fa file

infile=$1
rundir=`dirname $0`
binpath=$rundir

if [ "$infile" == "" ]; then 
    echo "usage:  sortedTopoMSA2numTMbardef.sh infile"
    exit
fi

red=#FF0000
green=#00FF00
blue=#0000FF
yellow=#FFFF00
pink=#FF00FF
cyan=#00FFFF
black=#000000


tmpfile=$(mktemp /tmp/tmp.sortedTopoMSA2numTMbardef.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  

python $binpath/countTMregion.py $infile | sort -k2,2rg > $tmpfile
maxNumTM=`head -n 1 $tmpfile | awk '{print $2}'`
minNumTM=`tail -n 1 $tmpfile | awk '{print $2}'`
awk -v maxNumTM=$maxNumTM -v minNumTM=$minNumTM '
{
sub("^>*", "");
seqid=$1;
if (seqid != "Consensus"){
    print seqid "," $2
}
}
END{
prin ""
}
' $tmpfile

rm -f $tmpfile
