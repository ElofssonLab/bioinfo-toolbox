#!/bin/bash

progname=`dirname $0`
usage="
Usage: $progname fasta-seq-file [-outpath outpath]
"
PrintHelp() {
    echo "$usage"
}
binpath=$DATADIR3/bin
infile=
outpath=./

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
            -q|-quiet|--quiet) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

nodename=`uname -n`
case $nodename in 
    ferlin*pdc.kth.se) 
        echo "OK, you are running from pdc node" 
        ;;
    *) 
        echo "Error. This script should be run from PDC login node. Exit" >&2
        exit 1
        ;;
esac

if [ "$infile" == ""  ]; then 
    echo "infile not set. Exit" >&2
    exit 1
elif [ ! -s "$infile" ] ; then
    echo "infile '$infile' not exist or empty. Exit. " >&2
else 
    /bin/mkdir -p $outpath
    /bin/mkdir -p $outpath/log
    echo "Command line: $commandline"
    echo "Start at    `/bin/date`"
    res1=$(/bin/date +%s.%N)
    basename=`basename $infile`
    rootname=${basename%.*}

    numseq=`$binpath/countseq.py $infile -nf `
    if [ "$numseq" == "" ]; then 
        echo "Fatal. countseq.py error. Exit."  >&2
        exit 1
    elif [ $numseq -le 0 ] ; then
        echo "Zero sequence found in file '$infile'. Exit."  >&2
        exit 1
    else
        echo "$numseq sequences are going to be predicted by topcons-single."
        numseqperjob=`expr 3 \* 3600 \* 2`
        cnt=0
        for ((i=0;i<=numseq;i+=numseqperjob)); do
            ((beg=i))
            ((end=i+numseqperjob))
            splittedfastafile=$outpath/${rootname}_split${cnt}.fa
            logfile=$outpath/log/${rootname}_split${cnt}.log
            #echo "$binpath/catfasta.py -b $beg -e $end $infile > $splittedfastafile"
            $binpath/catfasta.py -b $beg -e $end $infile > $splittedfastafile
            if [ -s "$splittedfastafile" ]; then 
                #echo "$splittedfastafile"
                esubmit -n 1 -t 240 "$DATADIR3/wk/MPTopo/bin/run_topcons_single.sh $splittedfastafile -outpath $outpath > $logfile 2>&1"
            fi
            ((cnt++))
        done
        echo "Finished at `/bin/date`"
    fi
fi
