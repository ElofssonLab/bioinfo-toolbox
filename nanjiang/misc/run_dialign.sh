#!/bin/bash

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname [-l LISTFILE] -outpath DIR  [-q]
        $wspace FILE 
Options:
  -outpath DIR   Set output path
  -l       FILE  Set the fileListFile, one filename per line
  -force, -overwrite 
                 Overwirte the result
  -maxseq INT    Maximal number of sequence to run, default: 1000
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2013-05-14, updated 2013-05-14, Nanjiang Shu 
"
dialign_exe=dialign2-2

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
Run_DIALIGN(){ #{{{
    local seqfile="$1"
    if [ ! -s "$seqfile" ]; then
        echo "seqfile $seqfile is empty. Ignore" >&2
        return 1
    fi

    local numseq=`grep "^>" "$seqfile" | wc -l`

    if [ $numseq -gt $maxseq ]; then 
        echo "numseq $numseq is over $maxseq. Ignore $seqfile"
        return 1
    fi

    local basename=`basename $seqfile`
    local rootname=${basename%.*}
    local tmpseqfile=$outpath/$basename.tmp
    local tmpalifile=$outpath/$basename.tmp.ali
    local tmpmfafile=$outpath/$basename.tmp.fa
    local alifile=$outpath/$rootname.dialign.ali
    local mfafile=$outpath/$rootname.dialign.mfa
    if [ -s "$mfafile" -a $isOverwrite -eq 0 ] ;then 
        echo "result $mfafile exist. Ignore $seqfile"
        return 1
    fi

    # run dialign
    /bin/cp -uf $seqfile $tmpseqfile
    $dialign_exe -fa $tmpseqfile 

    if [  -f "$tmpalifile" ] ; then
        mv $tmpalifile $alifile
    fi

    if [  -f "$tmpmfafile" ] ; then
        mv $tmpmfafile $mfafile
        echo $mfafile created
    fi

    rm -f $tmpseqfile

} 
#}}}


if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=0
outpath=
fileListFile=
fileList=()
maxseq=1000
isOverwrite=0

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
            -outpath|--outpath) outpath=$2;shift;;
            -l|--l|-list|--list) fileListFile=$2;shift;;
            -force|--force|-overwrite|--overwrite) isOverwrite=1;;
            -maxseq|--maxseq) maxseq=$2;shift;;
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

IsProgExist $dialign_exe

if [ "$outpath" == "" ]; then
    echo "outpath not set. exit" >&2
    exit 1
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi


for ((i=0;i<numFile;i++));do
    file=${fileList[$i]}
    Run_DIALIGN "$file"
done

