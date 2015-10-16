#!/bin/bash

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname [-l LISTFILE] [-outpath DIR]  [-q]
        $wspace FILE 
Options:
  -outpath DIR      Set output path
  -datapath DIR     Set datapath
  -maxdistkr INT    Set max distance for KR, default: 50
  -o       OUTFILE  Set output file
  -l       FILE     Set the idListFile, one filename per line
  -q                Quiet mode
  -h, --help        Print this help message and exit

Created 2009-05-26, updated 2012-09-19, Nanjiang Shu
Example: 
   $progname -l idlistfile -datapath data  -outpath outdir
"
PrintHelp(){ #{{{
    echo "$usage"
}
#}}}
rel2abs(){  #{{{ #fileordir
    #convert a path to absolute path, 
    #return "" if the file or path does not exist
    if [ ! -e "$1" ]; then 
        return 1
    fi
    curr="$PWD"
    if [ -d "$1" ]; then
        cd "$1" && echo "$PWD"       
    else 
        file="$(basename "$1")"
        cd "$(dirname "$1")" && dir="$PWD"
        echo "$dir/$file"
    fi
    cd "$curr"
    return 0
}
#}}}
AddAbsolutePath(){ #$fileordir#{{{
    #convert a path to absolute path, 
    #return "" if the file or path does not exist
    if [ ! -e "$1" ]; then 
        return 1
    fi
    curr="$PWD"
    if [ -d "$1" ]; then
        cd "$1" && echo "$PWD"       
    else 
        file="$(basename "$1")"
        cd "$(dirname "$1")" && dir="$PWD"
        echo "$dir/$file"
    fi
    cd "$curr"
    return 0
}
#}}}
IsProgExist(){ #{{{
    # usage: IsProgExist prog
    # prog can be both with or without absolute path
    type -P $1 &>/dev/null \
        || { echo The program \'$1\' is required but not installed. \
        Aborting $0 >&2; exit 1; }
    return 0
}
#}}}
IsPathExist(){ #{{{
# supply the effective path of the program 
    if ! test -d "$1"; then
        echo Directory \'$1\' does not exist. Aborting $0 >&2
        exit 1
    fi
}
#}}}
AnaKRBias(){ #{{{
    local id=$1
    aligned_seq_file=$datapath/$id.homology.cleaned.le1000.kalignP.fasta
    aligned_topo_file=$datapath/$id.sorted.orig.topomsa.fa
    seqfile=$outpath/$id.aaseq.fa
    topofile=$outpath/$id.topo
    sed 's/-//g' $aligned_topo_file > $topofile
    sed 's/-//g' $aligned_seq_file > $seqfile
    echo "/data3/wk/MPTopo/src/test_calculate_KR_bias.py -maxdist $maxDistKR -topofile $topofile -seqfile $seqfile"
    /data3/wk/MPTopo/src/test_calculate_KR_bias.py \
        -maxdist $maxDistKR \
        -topofile $topofile -seqfile $seqfile
} 
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=0
outpath=./
outfile=
idListFile=
idList=()

datapath=data
maxDistKR=50

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

if [ "$idListFile" != ""  ]; then 
    if [ -s "$idListFile" ]; then 
        while read line         
        do         
            idList+=("$line")
        done < $idListFile
    else
        echo listfile \'$idListFile\' does not exist or empty. >&2
    fi
fi

numFile=${#idList[@]}
if [ $numFile -eq 0  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi

if [ ! -d "$outpath" ] ; then
    mkdir -p $outpath
fi

for ((i=0;i<numFile;i++));do
    id=${idList[$i]}
    AnaKRBias "$id"
done

