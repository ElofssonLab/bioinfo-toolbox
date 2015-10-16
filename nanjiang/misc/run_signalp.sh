#!/bin/bash
# Filename: run_signalp.sh 

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname -t STR [-l LISTFILE] [-outpath DIR]  [-q]
        $wspace FILE 
Description:
    Run signalp prediction for the given sequence file
Options:
  -t       STR      Type of protein, gram-, gram+ or euk. (default: euk)
  -outpath DIR      Set output path
  -l       FILE     Set the fileListFile, one filename per line
  -q                Quiet mode
  -h, --help        Print this help message and exit

Created 2013-10-30, updated 2013-10-30, Nanjiang Shu
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
RunSignalp(){ #{{{
    local infile="$1"
    local basename=`basename "$infile"`
    local rootname=${basename%.*}
    local dirname=`dirname "$infile"`
    local outpath=
    if [ "$g_outpath" == "" ];then
        outpath=$dirname
    else
        outpath=$g_outpath
    fi
    local outfile=$outpath/$rootname.signalp_short
    exec_cmd "$exec_signalp -t $type -f short $infile > $outfile"
} 
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=0
g_outpath=
fileListFile=
fileList=()
exec_signalp=$DATADIR3/usr/share/signalp-4.0/signalp
type=euk

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        fileList+=("$1")
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) g_outpath=$2;shift;;
            -t|--t) type=$2;shift;;
            -l|--l|-list|--list) fileListFile=$2;shift;;
            -q|-quiet|--quiet) isQuiet=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        fileList+=("$1")
    fi
    shift
done

IsProgExist "$exec_signalp"

if [ "$fileListFile" != ""  ]; then 
    if [ -s "$fileListFile" ]; then 
        while read line
        do
            fileList+=("$line")
        done < $fileListFile
    else
        echo listfile \'$fileListFile\' does not exist or empty. >&2
    fi
fi

numFile=${#fileList[@]}
if [ $numFile -eq 0  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi

if [ "$g_outpath" != "" -a ! -d "$g_outpath" ];then
    mkdir -p "$g_outpath"
fi

for ((i=0;i<numFile;i++));do
    file=${fileList[$i]}
    RunSignalp "$file"
done
