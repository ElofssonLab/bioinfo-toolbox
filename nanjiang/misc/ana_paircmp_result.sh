#!/bin/bash

usage="
usage:  $0 paircmp-file [paircmp-file ...] 
           [-outpath DIR]

Description: analyze the result of pairwise topology comparison

    -outpath DIR   Set dir for output result, (default: ./result)

Created 2011-11-02, updated 2011-11-07, Nanjiang Shu
"
if [ "$1" == "" ] ; then 
    echo "$usage"
    exit
fi
commandline="$0 $*"

binpath=/data3/wk/MPTopo/bin
paircmpFileList=()
outpath=result

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        paircmpFileList+=("$1")
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) echo "$usage"; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -nc|-not-clean|--not-clean) isClean=0;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        paircmpFileList+=("$1")
    fi
    shift
done

numFile=${#paircmpFileList[@]}
if [ $numFile -lt 1 ]; then
    echo "paircmpFile not set. Exit." >&2
    exit 1;
fi

for ((i=0;i<numFile;i++)); do
    paircmpFile="${paircmpFileList[$i]}"
    if [ ! -s "$paircmpFile" ] ; then 
        echo "paircmpFile '$paircmpFile' does not exist or empty. Ignore." >&2
        continue
    fi

    basename=`basename "$paircmpFile"`
    rootname=${basename%.*}
    mkdir -p "$outpath"

    # 1. Generate SeqIDT vs CmpClass table
    seqidtCmpClassTable="$outpath/$rootname.seqidt-cmpclass.table.txt"
    cmpclass="OK SHIFT INV INV_SHIFT DIFF"
awk '
BEGIN{
    split("OK SHIFT INV INV_SHIFT DIFF", cmpclass, " ");
    for(i=0; i < 10; i++){
        for (j = 1; j <= 5;j++){
            freq[i, cmpclass[j]] = 0;
        }
        subsum[i] = 0;
    }
}
/Pairwise/{
    for (i=0;i <10;i++) {
        if($4>=i*10 && $4 <(i+1)*10) {
            freq[i,$5] ++;
            subsum[i] ++;
        }
    }
}
END{
    printf("%4s %s\t","#Idx", "SeqIDT");
    for (j=1;j<=5;j++){
        printf( " %10s", cmpclass[j]);
    }
    printf(" %10s","Occurrence");
    printf ("\n");
    for(i=0;i<10;i++){
        printf("%4d %d-%d\t",i, i*10,(i+1)*10);
        for (j=1;j<=5;j++){
            printf (" %10.3f", freq[i,cmpclass[j]]/subsum[i]*100);
        }
        printf(" %10d",subsum[i]);
        printf ("\n");
    }
}
' "$paircmpFile" > "$seqidtCmpClassTable"

    echo "Result output to $seqidtCmpClassTable"

    echo "Draw figure for $seqidtCmpClassTable"

    $binpath/plotCmpClass.sh "$seqidtCmpClassTable" -o eps  -outpath "$outpath"
done 
