#!/bin/bash 

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname -l fafilelist [FILE [FILE...]]
Options:
  -o       OUTFILE  Set output file
  -q                Quiet mode
  -h, --help        Print this help message and exit

Created 2012-05-28, updated 2012-05-28, Nanjiang Shu  
"
PrintHelp(){ #{{{
    echo "$usage"
}
#}}}

CalPercentTMPro(){ # $fafile
    local fafile="$1"
    if [ ! -s "$fafile"  ]; then 
        echo seqfile $fafile does not exist or empty. Ingore.
        return 1
    fi

    tmpfile=$(mktemp /tmp/tmp.calPercentTMPro_scampi.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }   
    trap 'rm -rf "$tmpfile"' INT TERM EXIT

    $DATADIR3/program/newscampiscript/isMemPro_scampi.sh "$fafile" > $tmpfile

    if [ ! -s $tmpfile ]; then 
        echo isMemPro_scampi.sh failed for file $fafile. Ignore
        rm -f $tmpfile
        return 1
    fi

    rst=`awk 'BEGIN{y=0;n=0;total=0}{if($2=="yes"){y++}else if($2=="no"){n++} total++}END{printf("%6d / %6d = %4.1f", y, total, y/total*100)}' $tmpfile`
    basename=`basename $fafile`
    rootname=${basename%.*}
    echo $rootname $rst
    rm -f $tmpfile
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=0
outfile=
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
            -o|--o) outfile="$2";shift;;
            -l|--l|-list|--list) fileListFile="$2";shift;;
            -q|-quiet|--quiet) isQuiet=1;;
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
    echo Input not set! Exit. >&2
    exit 1
fi

for ((i=0;i<numFile;i++));do
    file=${fileList[$i]}
    CalPercentTMPro "$file"
done

