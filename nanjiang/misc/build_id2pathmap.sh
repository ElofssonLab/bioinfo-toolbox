#!/bin/bash
progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname DIR [-o OUTFILE]
Description:
    build id to subdir path map, where real data is stored under subfolders
    split_0 split_1...
Options:
  -o   OUTFILE   Output the result to file
  -p   STR       search pattern
  -h             Print this help message and exit

Created 2013-06-14, updated 2013-06-14, Nanjiang Shu 

Examples:
  $progname hhprofile -p \"*.hhm\" > hhprofile/id2pathmap.txt
  $progname hhalign -p \"*.hhr\" > hhalign/id2pathmap.txt
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
MyFunc(){ #{{{
    echo "\"$1\""
} 
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

datapath=
outfile=
pattern=

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        datapath="$1"
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -o | --o) outfile="$2"; shift;;
            -p|--p|-pattern|--pattern) pattern="$2";shift;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        datapath="$1"
    fi
    shift
done

if [ "$outfile" == "" ];then
    outfile=/dev/stdout
fi

if [ "$datapath" == "" ];then
    echo  datapath not set, exit >&2
    exit
fi
if [ "$pattern" == "" ];then
    echo  pattern not set, exit >&2
    exit
fi

dirlist=`find "$datapath" -maxdepth 1 -type d | grep split_ | sort  -t"_" -k2,2g | rootname`
(for folder in $dirlist; do
    find $datapath/$folder -name "$pattern" | rootname |sort -u |awk -v FOLDER=$folder '{print $1 "\t" FOLDER}'
done) > "$outfile"


