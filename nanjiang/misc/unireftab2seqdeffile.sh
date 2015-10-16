#!/bin/bash

unireftabfile=$1

if [ "$unireftabfile" == "" ] ; then 
    echo "usage:  unireftab2seqdeffile.sh tabfile"
    exit 1
fi

awk -F"\t" '{
if (NR==1) {
print "#SeqID\tSeqName\tSeqLength\tOrganism";
}
else{
seqid=$1;
name=$3;
sub(/Cluster: /, "", name);
organism=$6;
len=$7;
print seqid "\t" name "\t" len "\t" organism;
}

}' $unireftabfile


