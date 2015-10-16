#!/bin/bash

psipreddir=$PSIPRED/../psipred32Fast
binpath=/data3/wk/MPTopo/bin

usage="
Usage:  run_kalignP.sh fasta-seq-file

run kalignP, add the posotion specific gap penalties estimated from
secondary structures predicted by PSIPRED_single

Options:
    -psipreddir <dir> : set the path for executables, default=$psipreddir 
    -q                : quite mode, do not write messages

    -h|--help         : print this help message and exit

Created 2010-12-03
nanjiang.shu@gmail.com

Examples: 
    run_kalignP.sh test.fasta
"
function IsProgExist()#{{{
# usage: IsProgExist prog
# prog can be both with or without absolute path
{
    type -P $1 &>/dev/null || { echo "The program \"$1\" is required but it's not installed. Aborting." >&2; exit 1; }
}
#}}}
function PrintHelp()
{
    echo "$usage"
}
function RunKalignP() #$file#{{{
{
    local file=$1
    local basename=`basename $file`
    local rootname=${basename%.*}
    local tmpdir=/tmp/$rootname$$
    mkdir -p $tmpdir
    $binpath/splitfasta.py $file -nameseq -outpath $tmpdir $quietOptionStr 2> $errFile
    if [ -s $errFile ]; then 
        cat $errFile  >&2                                          
        exit 1;
    fi

    local aaFileListFile=$tmpdir/aafilelist.txt
    find $tmpdir -name "*.aa" > $aaFileListFile

    #predict secondary structure
    
    $psipreddir/runpsipred_single $quietOptionStr -l $aaFileListFile -outpath $tmpdir/$rootname

    #add gap penalties
    $binpath/addGapPenaltyByPredSS.py $quietOptionStr -p-shift-c $p_shift_C -p-shift-he $p_shift_HE -p-threshold-c $p_threshold_C -p-threshold-he $p_threshold_HE -weight-c $weight_C -weight-he $weight_HE $file -sspath $tmpdir -outpath $tmpdir 2> $errFile
    if [ -s $errFile ]; then 
        cat $errFile >&2
        exit 1;
    fi
    gpFastaFile=$tmpdir/$rootname.gp.fa
    alnFile=$outpath/$rootname.kalignP.$format
    $binpath/kalignP $gpFastaFile $quietOptionStr -format $format -o $alnFile
    if [ "$isQuiet"  != "true" ]; then
        echo "The alignment has been output to $alnFile"
    fi

    rm -rf $tmpdir
}
#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

fileList=
listFile=
errFile=$(mktemp /tmp/tmperr.XXXXXXX) || { echo "Failed to create temp file"; exit 1; } 
outpath=./
isQuiet=false
quietOptionStr=
format=msf

p_shift_C=0.50
p_shift_HE=0.53
p_threshold_C=0.0
p_threshold_HE=0.0
weight_C=0.60
weight_HE=0.38

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -psipreddir|--psipreddir) psipreddir=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -l|--list) listFile=$2;shift;;
            -f|-foramt|--format) foramt=$2;shift;;
            -q|--q|-quiet|--quiet) isQuiet=true;quietOptionStr=-q;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        fileList="$fileList $1 "
    fi
    shift
done

if [ "$fileList"  == "" -a "$listFile" == "" ]; then 
    echo "Error! not input is set. Exit..."  >&2
    exit
fi
if [ ! -d "$psipreddir" ]; then 
    echo "Error!  psipreddir= \"$psipreddir\" does not exist. Exit..." >&2
    exit
fi

mkdir -p $outpath

if [ "$fileList" != "" ]; then 
    for file in $fileList; do 
        RunKalignP $file
    done
fi

if [ -s "$listFile" ]; then 
    for file in $(cat $listFile); do 
        RunKalignP $file
    done
fi 

rm -f $errFile
