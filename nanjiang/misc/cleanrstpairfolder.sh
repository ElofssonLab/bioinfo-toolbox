#!/bin/bash 
progname=`basename $0`
usage="
usage:  $progname DIR [DIR ...]
Description:
    Clean rstpair pairwise topology comparison folder

OPTIONS:
  -h, --help     Print this help message and exit

Created 2013-05-08, updated 2013-05-08, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

CleanRSTPAIRFolder(){ ##{{{
    local folder="$1"
    local tmpfile=
    #trap 'rm -f "$tmpfile"' INT TERM EXIT

    find "$folder" -type f -name "*.topoaln.fa" -print0 | xargs -0 rm -f
    find "$folder" -type f -name "*.eps" -print0 | xargs -0 rm -f
    find "$folder" -type f -name "*.topoaln.topowithdgscore" -print0 | xargs -0r gzip -Nf
    find "$folder" -type f -name "*.dgscorelist" -print0 | xargs -0r gzip -Nf
    find "$folder" -type f -name "*.paircmp" -print0 | xargs -0r gzip -Nf
    find "$folder" -type f -name "*.txt" -print0 | xargs -0r gzip -Nf
    echo "$folder cleaned"
}
#}}}

folderList=()

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        folderList+=("$1")
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        folderList+=("$1")
    fi
    shift
done

numFolder=${#folderList[@]}
if [ $numFolder -eq 0  ]; then
    echo Input not set! Exit. >&2
    PrintHelp
    exit 1
fi

for ((i=0;i<numFolder;i++));do
    folder=${folderList[$i]}
    CleanRSTPAIRFolder "$folder"
done
