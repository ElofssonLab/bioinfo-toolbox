#!/bin/bash

#binpath=/data3/wk/MPTopo/bin
binpath=$BINPATH
embossbin=/usr/bin

if [ "$EMBOSSBIN" != "" ]; then
    embossbin=$EMBOSSBIN
fi

usage="
Usage:   calEvoDist.sh fasta-seq-file
Compute the evolutionary distance of two sequences

Options:
    -l <file> : list file for fasta sequence files
    -o <file> : write the output to file
    -f 0|1    : 0 for raw seq, 1 for aligned seq, default=0
    -emboss
    -q        : quite mode, do not write messages
    -h|--help : print this help message and exit

Created 2011-04-29, updated 2012-09-27, Nanjiang Shu  
nanjiang.shu@gmail.com

Examples: 
    calEvoDist.sh test.fasta
"
PrintHelp() {
    echo "$usage"
}
IsProgExist(){ #{{{
# usage: IsProgExist prog
# prog can be both with or without absolute path
    type -P $1 &>/dev/null || { echo "The program \"$1\" is required but it's not installed. Aborting." >&2; exit 1; }
}
#}}}
CalEvoDist(){ #$file#{{{
    local file=$1
    local basename=`basename $file`
    local rootname=${basename%.*}
    local tmpdir=$(mktemp -d /tmp/tmpdir.calEvoDist.XXXXXXX) || { echo "Can not create tmpdir" >&2; exit; }
    local alnfile=$tmpdir/$rootname.aln.fasta
    if [ $format -eq 0 ]; then 
        $binpath/kalign -q -f fasta $file -o $alnfile 2> /dev/null
    else
        /bin/cp -f $file $alnfile
    fi

    # duplicate sequence for tree-puzzle, since tree-puzzle require >=3
    # sequences
    local alnfile_dup=$tmpdir/$rootname.alndup.fasta
    /bin/cat $alnfile $alnfile | /usr/bin/awk 'BEGIN{N=0}{if (substr($0,1,1)==">"){N++; if (N==3){print ">dupseq1"}else if (N==4){print ">dupseq2"}else {print} } else {print}}' > $alnfile_dup
    # convert the sequence format to PHYLIP
    local alnfile_phylip=$tmpdir/$rootname.aln.phy
    $embossbin/seqret -sf fasta -sequence $alnfile_dup -osf phylip -outseq $alnfile_phylip >/dev/null 2>&1

    # calculate distance
    echo "y" > $tmpdir/yes
    $binpath/treepuzzle  $alnfile_phylip $tmpdir/$rootname > /dev/null 2>&1 < $tmpdir/yes
    distfile=$tmpdir/$rootname.dist
    #cat $distfile                
    /usr/bin/awk 'BEGIN{id1=""; id2=""; dist=0.0}{if(NR==2){id1=$1;dist=$3} if(NR==3){id2=$1}} END{print id1,id2,dist}' $distfile >> $outfile

    rm -rf $tmpdir
}
#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

fileList=()
fileListFile=
outpath=./
isQuiet=false
quietOptionStr=
format=0
outfile=

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        fileList+=("$1")
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -psipreddir|--psipreddir) psipreddir=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -o|-outfile|--outfile) outfile=$2;shift;;
            -l|--list) fileListFile=$2;shift;;
            -f|-foramt|--format) foramt=$2;
                if [ "$2" != "0" -a "$2" != "1" ]; then
                    echo "Error! format can only be 0 or 1. Aborting." >&2; 
                    exit
                fi
                shift;;
            -q|--q|-quiet|--quiet) isQuiet=true;quietOptionStr=-q;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        fileList+=("$1")
    fi
    shift
done

IsProgExist $binpath/kalign
IsProgExist $binpath/treepuzzle
IsProgExist $embossbin/seqret

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

if [ "$outfile" == "" ]; then 
    outfile=/dev/stdout
else
    cat /dev/null > $outfile
fi

for ((i=0;i<numFile;i++));do
    file=${fileList[$i]}
    CalEvoDist $file
done

