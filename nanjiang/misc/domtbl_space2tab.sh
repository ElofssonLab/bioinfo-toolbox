#!/bin/bash
# convert space to tab for domtbl files to reduce the file size

usage="
USAGE: $0 domtblfile [domtblfile ...] [ -l LISTFILE ]

Created 2014-12-12, updated 2014-12-12, Nanjiang Shu
"
exec_cmd(){
    echo "$*"
    eval "$*"
}

PrintHelp(){ #{{{
    echo "$usage"
}
#}}}

Domtbl_space2tab(){
    local infile=$1
    local tmpfile=${infile}.bak.$$
    sed -e  '/^[^#]/ s/\s\+/\t/g' -e 's/[ \t]*$//' $infile > $tmpfile
    /bin/mv -f $tmpfile $infile
    echo "$infile space converted to tab"
}

if [ $# -lt 1 ];then
    echo "$usage"
    exit 1
fi
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

fileListFile=
fileList=()

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
            -l|--l|-list|--list) fileListFile=$2;shift;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        fileList+=("$1")
    fi
    shift
done

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
    echo $0: Input not set! Exit. >&2
    exit 1
fi

for ((i=0;i<numFile;i++));do
    file=${fileList[$i]}
    Domtbl_space2tab "$file"
done

