#!/bin/bash

idlistfile=idlist.txt
datapath=data
outpath=out2
maxDistKR=50

idList=()

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        idList+=("$1")
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -datapath|--datapath) datapath=$2;shift;;
            -maxdistkr|--maxdistkr) maxDistKR=$2;shift;;
            -o|--o) outfile=$2;shift;;
            -l|--l|-list|--list) idListFile=$2;shift;;
            -q|-quiet|--quiet) isQuiet=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        idList+=("$1")
    fi
    shift
done

if [ ! -d $outpath ]; then
    mkdir -o $outpath
fi


for id in $(cat $idlistfile); do
    echo "======================================"
    echo "$id"
    echo "
/data3/wk/MPTopo/src/drawMSATopo.py \
$datapath/$id.sorted.orig.topomsa.fa \
-aaseq $datapath/$id.homology.cleaned.le1000.kalignP.fasta \
-krbias -outpath $outpath
    "
    /data3/wk/MPTopo/src/drawMSATopo.py \
        $datapath/$id.sorted.orig.topomsa.fa \
        -aaseq $datapath/$id.homology.cleaned.le1000.kalignP.fasta \
        -krbias -outpath $outpath -maxdistkr $maxDistKR
    figfile=$outpath/$id.sorted.orig.topomsa.png
    thumbnail=$outpath/thumb.$id.sorted.orig.topomsa.png
    if [ -f $figfile ]; then
        convert -thumbnail 200 $figfile $thumbnail
    fi
done
