#!/bin/bash

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname DIR [DIR...]
Description:
    rename result file for ana_paircmp_result.py
Options:
  -h, --help        Print this help message and exit

Created 2013-07-04, updated 2013-07-04, Nanjiang Shu 
"
PrintHelp(){ #{{{
    echo "$usage"
}
#}}}

RenameFileName(){
    local dir="$1"
    currdir=$PWD
    cd $dir
    for file in *type_all_dg100.0_gap0.0_ps0.0*; do newname=`echo $file| sed 's/type_all_dg100.0_gap0.0_ps0.0//g'`; mv $file $newname; done
    cd $currdir
}


if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

dirList=()

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        dirList+=("$1")
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        dirList+=("$1")
    fi
    shift
done


numDir=${#dirList[@]}
if [ $numDir -eq 0  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi

for ((i=0;i<numDir;i++));do
    dir=${dirList[$i]}
    RenameFileName "$dir"
done

