#!/bin/bash

psipreddir=$PSIPRED/../psipred32
usage="
Usage:  my_runpsipred.sh fasta-seq-file

run psipred prediction (psiblast will be invokded to get the pssm)
multiple sequences can be included in the fasta-seq-file
the predicted secondary structure data files will be output
to \$outpath/\$rootname(fasta-seq-file) and named as \$rootname_i, i=0,1,2...

Options:
    -l          <file>: set the list file for fasta-seq-file
    -psipreddir <dir> : set the path for executables, default=$psipreddir 
    -outpath    <dir> : output the result to dir, default = ./
    -q                : quite mode, do not write messages
    -h|--help         : print this help message and exit

Created 2010-10-22, updated 2010-10-22
nanjiang.shu@gmail.com

Examples: 
    my_runpsipred.sh test.fasta
"
function PrintHelp()
{
    echo "$usage"
}
function MyRunPsiPred() #$file#{{{
{
    local file=$1
    local basename=`basename $file`
    local rootname=${basename%.*}
    local suboutpath=$outpath/$rootname
    mkdir -p $suboutpath
    splitfasta.py $file -nameseq -outpath $suboutpath -q 2> $errFile
    if [ -s $errFile ]; then 
        cat $errFile  >&2
        exit 1;
    fi
    currDir=$PWD
    cd $suboutpath
    local aaFileListFile=aafilelist.txt
    find ./ -name "*.aa" > $aaFileListFile
    (for aaFile in $(cat $aaFileListFile); do 
         $psipreddir/runpsipred $aaFile 
      echo
    done) 1> /dev/null 2> $errFile
    if [ -s $errFile ]; then 
        cat $errFile >&2
        exit 1;
    fi
    numseq=`wc -l $aaFileListFile | awk '{print $1}'`
    if [ "$isQuiet" != "true" ] ; then 
        echo -e "$numseq \tsequences in the file $file have been predicted by psipred"
    fi
    # merge the separated prediction into one
    mergedSS2File=${rootname}_all.ss2
    (for ((i=0;i<$numseq;i++));do
        ss2File=${rootname}_$i.ss2
        echo "#PSIPRED $ss2File"
        cat $ss2File
    done ) > $mergedSS2File
    cd $currDir
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
            -q|--q|--quiet) isQuiet=true;;
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
        MyRunPsiPred $file
    done
fi

if [ -s "$listFile" ]; then 
    for file in $(cat $listFile); do 
        MyRunPsiPred $file
    done
fi 

rm -f $errFile
