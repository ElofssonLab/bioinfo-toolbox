#!/bin/bash
progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname HHprofileDir
Options:
  -h, --help        Print this help message and exit

Created 2012-12-03, updated 2012-12-03, Nanjiang Shu 
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

hhprofileDIR=

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        hhprofileDIR=$1
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        hhprofileDIR=$1
    fi
    shift
done

if [ "$hhprofileDIR" == "" ];then
    echo  hhprofileDIR not set, exit >&2
    exit
fi

dirlist=`find $hhprofileDIR -maxdepth 1 -type d | grep split_ | sort  -t"_" -k2,2g | rootname`
for folder in $dirlist; do
    find $hhprofileDIR/$folder -name "*.hhm" | rootname |sort -u |awk -v FOLDER=$folder '{print $1 "\t" FOLDER}'
done


