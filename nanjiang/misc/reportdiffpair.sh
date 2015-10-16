#!/bin/bash

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname diffpairfile [-outpath DIR]  [-q]
        $wspace FILE 
Options:
  -outpath DIR      Set output path
  -seqdb   FILE     Set sequence db file
  -topodb  FILE     Set topology db file
  -o       OUTFILE  Set output file
  -l       FILE     Set the fileListFile
  -q                Quiet mode
  -h, --help        Print this help message and exit

Created 2009-05-26, updated 2011-11-19, Nanjiang Shu
"
PrintHelp(){ #{{{
    echo "$usage"
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
ReportDiffPair() {
    local diffpairfile=$1
    tmpdir=$(mktemp -d /tmp/tmpdir.reportdiffpair.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }  
    pairlistfile=$tmpdir/t1.pairlist
    awk '{}' $diffpairfile > $pairlistfile


}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=0
outpath=./
outfile=
diffpairfile=

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        diffpairfile=$1
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -i | --i) diffpairfile=$2; shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -o|--o) outfile=$2;shift;;
            -q|-quiet|--quiet) isQuiet=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        diffpairfile=$1
    fi
    shift
done


if [ "$diffpairfile" == ""  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi

ReportDiffPair $diffpairfile
