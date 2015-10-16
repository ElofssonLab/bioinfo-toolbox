#!/bin/bash
# Description: get pair list from dumped tableinfo files.
# For each group, at most 10000 pairs are selected

usage="
usage dump_select_pairlist.sh tableinfo-file [tableinfo-file ...]
  -seqidttype INT    Set sequence identity type, (default: 0)
                     0: seqIDT = numIDTRes /alnLength
                     1: seqIDT = numIDTRes / min(len1, len2)
                     2: seqIDT = numIDTRes / (alnLength - NumGAP)

Created 2011-11-01, updated 2012-05-11, Nanjiang Shu  

Example:
    dump_select_pairlist.sh tableinfofiles ... -seqidttype 1 > outfile
"
filelist=
N_MAX=30000
seqidttype=0


isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        filelist="$filelist $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) echo "$usage"; exit;;
            -seqidttype|--seqidttype) seqidttype=$2;shift;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        filelist="$filelist $1"
    fi
    shift
done


if [ "$filelist" == "" ]; then 
    echo "$usage"
    exit 1;
fi


tmpfile=$(mktemp /tmp/tmp.dump_select_pairlist.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
dumpedfile=$(mktemp /tmp/tmp.dump_select_pairlist.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  

cat $filelist > $dumpedfile

SEQIDTIDX=3

case $seqidttype in 
    0) SEQIDTIDX=3 
        ;;
    1) SEQIDTIDX=12 
        ;;
    2) SEQIDTIDX=13 
        ;;
    *)
        echo Wrong seqidttype. Exit
esac

for ((i=0;i<10;i++)); do
    awk -v I=$i -v IDX=$SEQIDTIDX '{if($IDX>=I*10 && $IDX<(I+1)*10 && NF>10) if($1>$2){print $1, $2} else{print $2,$1}}'  $dumpedfile | sort -u | randlist | head -n $N_MAX >$tmpfile
    wc -l $tmpfile >&2
    cat $tmpfile
done

rm -f $dumpedfile
rm -f $tmpfile
