#!/bin/bash

# Created 2012-02-23, updated 2012-02-23, Nanjiang Shu  
# create color strips dataset for ITOL phylogenetic tree drawing from 
# sorted.orig.topomsa.fa file

infile=$1

if [ "$infile" == "" ]; then 
    echo "usage:  sortedTopoMSA2colordef.sh infile"
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
seqid=$1;
tag=$2;
if (tag=="IDT"){
    print seqid "," green
}
else if (tag == "OK"){
    print seqid "," red
}
else if (tag == "SHIFT"){
    print seqid "," pink
}
else if (tag == "INV"){
    print seqid "," blue
}
else if (tag == "INV_SHIFT"){
    print seqid "," cyan
}
else if (tag == "DIFF"){
    print seqid "," black
}
}
END{
print ""
}
' $infile
