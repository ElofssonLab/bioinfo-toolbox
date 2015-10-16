#!/bin/bash

# Created 2012-03-16, updated 2012-03-16, Nanjiang Shu 
# create color strips dataset for ITOL phylogenetic tree drawing from 
# clustered.orig.topomsa.fa

infile=$1

if [ "$infile" == "" ]; then 
    echo "usage:   clusteredTopoMSA2colordef.sh infile"
    exit
fi

red=#FF0000
green=#00FF00
blue=#0000FF
yellow=#FFFF00
pink=#FF00FF
cyan=#00FFFF
black=#000000


awk -v red=$red -v green=$green -v blue=$blue -v yellow=$yellow -v pink=$pink -v cyan=$cyan -v black=$black '
/^>/ {
sub("^>*", "");
seqid=$1
gsub(",","",seqid)
tag=$3;
if (tag=="ClusterNo=1"){
    print seqid "," "#0000FF"
}
else if (tag=="ClusterNo=2"){
    print seqid "," "#6363FF"
}
else if (tag=="ClusterNo=3"){
    print seqid "," "#8080FF"
}
else if (tag=="ClusterNo=4"){
    print seqid "," "#ABABFF"
}
else if (tag=="ClusterNo=5"){
    print seqid "," "#9B84FF"
}
else if (tag=="ClusterNo=6"){
    print seqid "," "#A0DFFF"
}
else if (tag=="ClusterNo=7"){
    print seqid "," "#DEECFF"
}
else{
    print seqid "," "#EDFFFF"
}
}
END{
print ""
}
' $infile
