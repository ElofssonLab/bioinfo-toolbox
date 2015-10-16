#!/bin/bash
# run pfam_scan.pl
# Note, the copy of data on clusters can be performed in the master script
# this script should keep clean
# specify pfamdb to local if running in clusters

step=1
nodename=`uname -n`
datadir3=$DATADIR3


pfamdb=$DATADIR/pfam/pfam27.0
path_pfamscan=$datadir3/usr/share/PfamScan/PfamScan
export PERL5LIB=${path_pfamscan}:$PERL5LIB
binpath=$datadir3/bin

which hmmscan

usage="
USAGE: $0 fasta-seq-file [-outpath DIR]
DESCRIPTION: 
    run pfam_scan.pl
    output file will be \$outpath/\$rootname.pfamscan

OPTIONS:
  -outpath DIR      Set output path, (default: ./)
  -pfamdb  STR      Set pfamdb, (default: $pfamdb)
  -r STR            Set options for pfam_scan.pl
                    e.g. -pfamb
  -q                Quiet mode
  -h,--help         Print this help message and exit

Created 2014-11-18, updated 2014-11-18, Nanjiang Shu
"
PrintHelp() {
    echo "$usage"
}
exec_cmd(){
    echo "$*"
    eval "$*"
}
function IsProgExist()#{{{
# usage: IsProgExist prog
# prog can be both with or without absolute path
{
    type -P $1 &>/dev/null || { echo "The program \"$1\" is required but it's not installed. Aborting $0" >&2; exit 1; }
}
#}}}
function IsPathExist()#{{{
# supply the effective path of the program 
{
    if ! test -d $1; then
        echo "Directory $1 does not exist. Aborting $0" >&2
        exit
    fi
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit 1
fi

commandline="$0 $*"

isQuiet=0
outpath=./
infile=
pfamscanOpt=

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        infile=$1
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -pfamdb|--pfamdb) pfamdb=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -r|--r) pfamscanOpt="$2"; shift;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

IsProgExist $binpath/countseq.py
IsProgExist $binpath/catfasta.py
IsPathExist $path_pfamscan


if [ "$infile" == ""  ]; then 
    echo "infile not set. Exit" >&2
    exit 1
elif [ ! -s "$infile" ] ; then
    echo "infile '$infile' not exist or empty. Exit. " >&2
else 
    mkdir -p $outpath
    mkdir -p $outpath/log
    echo "Command line: $commandline"
    echo "Start at    `/bin/date`"
    res1=$(/bin/date +%s.%N)
    basename=`basename $infile`
    rootname=${basename%.*}

    outfile=$outpath/$rootname.pfamscan
    outfile=`readlink -f $outfile`

    if [ -f $outfile ]; then 
        rm -rf $outfile
    fi
    numseq=`$binpath/countseq.py $infile -nf `
    if [ "$numseq" == "" ]; then 
        echo "Fatal. countseq.py error. Exit."  >&2
        exit 1
    elif [ $numseq -le 0 ] ; then
        echo "Zero sequence found in file '$infile'. Exit."  >&2
        exit 1
    else
        echo "$numseq sequences are going to be predicted by pfam_scan.pl."
        runningtime=`echo "$numseq / 2" | bc `
        echo "It takes about $runningtime seconds to run the prediction."
        exec_cmd "$path_pfamscan/pfam_scan.pl -fasta $infile -dir $pfamdb  -outfile $outfile $pfamscanOpt"
        echo "Finished at `/bin/date`"
        res2=$(/bin/date +%s.%N)
        printf "Running time for predicting %d sequences is: %.3F seconds\n" $numseq $(echo "$res2 - $res1"|/usr/bin/bc )
        echo "Result output to $outfile"
    fi
fi

