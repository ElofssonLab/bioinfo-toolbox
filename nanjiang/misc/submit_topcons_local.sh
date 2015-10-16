#!/bin/bash

run_tp_exe=$DATADIR3/wk/MPTopo/src/run_topcons.sh
progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname fafile  -outpath DIR
Options:
  -outpath DIR      Set output path
  -npara   INT      Number of parallel jobs, (default: NCPU/2)
  -q                Quiet mode
  -h, --help        Print this help message and exit

Created 2014-10-08, updated 2014-10-08, Nanjiang Shu 
"
PrintHelp(){ #{{{
    echo "$usage"
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
exec_cmd(){ #{{{
    echo "$*"
    eval "$*"
}
#}}}
RunTopcons(){ #{{{
    local fafile=$1
    local basename=`basename $fafile`
    local rootname=${basename%.*}
    local splittedpath=$outpath/splitted
    mkdir -p $splittedpath
    splitfasta.py $fafile -nfile $cpu -outpath $outpath/splitted -ext fa
    for ((i=0;i<cpu;i++)); do
        local splittedfafile=$outpath/splitted/${rootname}_$i.fa
        local splittedlogfile=$outpath/splitted/${rootname}_$i.log
        exec_cmd "$run_tp_exe $splittedfafile -outpath $splittedpath > $splittedlogfile 2>&1 &"
    done
    wait
    # dump the result
    (for ((i=0;i<cpu;i++)); do
        cat $outpath/splitted/${rootname}_$i.topcons.result.txt
    done) > $outpath/${rootname}.allinfo 
    echo Result output to  $outpath/${rootname}.topcons.result.txt
} 
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

numCPU=`cat /proc/cpuinfo | grep processor | wc -l`
if [ "$numCPU" != "" ];then
    cpu=`expr $numCPU / 2`
else
    cpu=1
fi
isQuiet=0
outpath=./
orgtype=

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        fafile=$1
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -npara|--npara) cpu=$2;shift;;
            -t|--t|-orgtype) orgtype=$2;shift;;
            -q|-quiet|--quiet) isQuiet=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        fafile=$1
    fi
    shift
done

if [ ! -f "$fafile" ]; then 
    echo fafile \'$fafile\' does not exist or empty. >&2
    exit
fi


if [ ! -d $outpath ] ; then 
    mkdir -p $outpath
fi
RunTopcons $fafile

