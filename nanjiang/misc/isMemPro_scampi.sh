#!/bin/bash

usage="
Usage:   isMemPro_scampi.sh fasta-seq-file
  check whether proteins are TM protein or not using scampi
Options:
  -o <file> : set outfile
  -q        : quiet mode
  -h|--help : print this help message and exit
Created 2011-03-08, updated 2011-03-08, Nanjiang
"
function PrintHelp()
{
    echo "$usage"
}
function IsProgExist()#{{{
# supply the effective path of the program 
{
    if ! test -x $1; then
        echo "program $1 does not exist, exit."
        exit
    fi
}
#}}}
function IsMemPro_Scampi()
{
    local fastafile=$1
    local tmpdir=`mktemp -d` 
    $binpath/mySCAMPI_run.pl -q $fastafile -outpath $tmpdir >/dev/null 2>&1
    basename=`basename $fastafile`
    xmlfile=$tmpdir/$basename.xml.res
    /usr/bin/awk -F: '{if($1=="SeqID"){printf ("%s\t", $2)}else if ($1=="IsTMProtein"){print $2}}' $xmlfile > $outfile
    /bin/rm -rf $tmpdir
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

binpath=$BINPATH
isQuiet=false
outfile=/dev/stdout
fastafile=

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        fastafile=$1
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -o|--o) outfile=$2;shift;;
            -q) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        fastafile=$1
    fi
    shift
done
IsProgExist $binpath/mySCAMPI_run.pl
IsProgExist /usr/bin/awk
IsProgExist /bin/rm

if [ "$fastafile" == "" ]; then
    echo "$0: Error, input not set!"
    exit
fi

if [ "$fastafile" != ""  ]; then
    if [ -f "$fastafile" ] ; then
        IsMemPro_Scampi $fastafile
    else 
        echo "seqfile $fastafile does not exist"
    fi
fi

