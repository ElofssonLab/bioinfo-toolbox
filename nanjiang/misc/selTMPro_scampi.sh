#!/bin/bash
usage="
Usage:   selTMPro_scampi.sh [-i] faFileList
  test argument passing
Options:
  -o       <file> : set output file
  -i       <file> : set the faFileList
  -th      <int>  : set threshold, default=75%
  -q              : quiet mode
  -h|--help       : print this help message and exit
Created 2011-08-15, updated 2011-08-15, Nanjiang
"
function PrintHelp()
{
    echo "$usage"
}
function AddAbsolutePath() #$path#{{{
{
    local var=$1
    if [ "${var:0:1}" != "/" ];then
        var=$PWD/$var # add the absolut path
    fi
    echo $var
    return 0
}
#}}}
function IsProgExist()#{{{
# usage: IsProgExist prog
# prog can be both with or without absolute path
{
    type -P $1 &>/dev/null || { echo "The program \"$1\" is required but it's not installed. Aborting." >&2; exit 1; }
}
#}}}
function IsPathExist()#{{{
# supply the effective path of the program 
{
    if ! test -d $1; then
        echo "Directory $1 does not exist. Aborting." >&2
        exit
    fi
}
#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outfile=/dev/stdout
faFileList=
threshold=75

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        fileList="$fileList $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -o|--outfile|-outfile) outfile=$2;shift;;
            -i) faFileList=$2;shift;;
            -th|--th) threshold=$2;shift;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        faFileList=$1
    fi
    shift
done

if [ "$faFileList" == "" ]; then
    echo "$0: Error, input not set!" >&2
    exit
fi

IsProgExist isMemPro_scampi.sh

for file in $(cat $faFileList); do 
    per=`isMemPro_scampi.sh $file | awk 'BEGIN{y=0;n=0;total=0}{if($2=="yes"){y++}else if($2=="no"){n++} total++}END{printf("%.0f" ,y/total*100)}'`
#     echo $per
#     echo $threshold
    if [ $per -ge $threshold ]; then
        echo "$file $per"
    fi
done
