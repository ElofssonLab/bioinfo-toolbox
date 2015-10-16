#!/bin/bash

#maxseq=1000
step=1
nodename=`uname -n`
TMPDIR=/scratch/$$
case $nodename in 
    *pdc.kth.se) 
    datadir3=$DATADIR3 
    TMPDIR=/scratch/$$
    ;;
    *uppmax.uu.se)
    datadir3=$DATADIR3
    if [ $SNIC_TMP ]; then
        TMPDIR=$SNIC_TMP
    else
        TMPDIR=/tmp/$$
    fi
    ;;
    *) 
    datadir3=/data3 
    if [ -d /scratch ]; then 
        TMPDIR=/scratch/$$
    else
        TMPDIR=/tmp/$$
    fi
    ;;
esac
if [ ! -d $TMPDIR ] ; then 
    mkdir -p $TMPDIR
fi


path_polyphobius=$datadir3/share/polyphobius
binpath=$datadir3/bin

# blastdb for polyphobius is created by
# $ ../blastget -ix uniref90.mem.fasta.ix -create uniref90.mem.fasta
# $ formatdb -i uniref90.mem.fasta -p T -o F  

usage="
usage:   run_polyphibius.sh fasta-seq-file [-outpath DIR]
Description: Run polyphibius
             output file will be $outpath/$rootname.polyphobius
  -outpath DIR      Set output path, (default: ./)
  -dbname STR       Set the blast dbname, default, uniref90.mem.fasta
  -q                Quiet mode
  -h|--help         Print this help message and exit
Created 2012-07-27, updated 2012-08-16, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}
IsProgExist(){ #{{{
# usage: IsProgExist prog
# prog can be both with or without absolute path
    type -P $1 &>/dev/null || { echo "The program \"$1\" is required but it's not installed. Aborting $0" >&2; exit 1; }
}
#}}}
IsPathExist(){ #{{{
# supply the effective path of the program 
    if ! test -d $1; then
        echo "Directory $1 does not exist. Aborting $0" >&2
        exit
    fi
}
#}}}
RunPolyPhobius(){ # $outdir#{{{
    local fafile=$1
    local outfile=$2
    local suboutdir=$3
    $path_polyphobius/run_phobius.sh $fafile $suboutdir $BLASTDB/$dbname $BLASTDB/$dbname.ix
    local bsname=`basename $fafile`
    local rtname=${bsname%.*}
    local rstfile=$suboutdir/$rtname.polyphobius
    echo rstfile $rstfile output
    cat $rstfile
    if [ -s $rstfile ]; then
        cat $rstfile >> $outfile
    fi
    return 0
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
dbname=uniref90.mem.fasta

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
            -dbname|--dbname) dbname=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done


IsProgExist $binpath/countseq.py
IsProgExist $binpath/catfasta.py
IsProgExist $binpath/splitfasta.py
IsPathExist $path_polyphobius

if [ ! -s "$path_polyphobius/data/$dbname.phr" ]; then 
    echo "blastdb $path_polyphobius/$dbname does not exist. Exit" >&2
    exit 1
fi

/bin/mkdir -p $outpath

if [ "$infile" == ""  ]; then 
    echo "infile not set. Exit" >&2
    exit 1
elif [ ! -s "$infile" ] ; then
    echo "infile '$infile' not exist or empty. Exit. " >&2
else 
    echo "Command line: $commandline"
    echo "Start at    `/bin/date`"
    res1=$(/bin/date +%s.%N)
    basename=`basename $infile`
    rootname=${basename%.*}

    outfile=$outpath/$rootname.polyphobius_long
    errfile=$outpath/$rootname.polyphobius.err
    logfile=$outpath/$rootname.polyphobius.log

    if [ -s $outfile ]; then 
        cat  /dev/null > $outfile
    fi
    if [ -s $logfile ]; then 
        cat  /dev/null > $logfile
    fi
    if [ -s $errfile ]; then 
        cat  /dev/null > $errfile
    fi
    numseq=`$binpath/countseq.py $infile -nf `
    if [ "$numseq" == "" ]; then 
        echo "Fatal. countseq.py error. Exit."  >&2
        exit 1
    elif [ $numseq -le 0 ] ; then
        echo "Zero sequence found in file '$infile'. Exit."  >&2
        exit 1
    else
        echo "$numseq sequences are going to be predicted by polyphobius."
        runningtime=`expr $numseq / 2 \* 60`
        echo "It takes about $runningtime seconds to run the prediction."
        case $nodename in 
            *.pdc.kth.se|*uppmax.uu.se)
                export BLASTDB=$TMPDIR/blastdb
                mkdir -p $BLASTDB
                echo "cp -f $path_polyphobius/data/$dbname.* $BLASTDB/"
                cp -f $path_polyphobius/data/$dbname.* $BLASTDB/
                ;;
            *)
                export BLASTDB=$path_polyphobius/data
                ;;
        esac
        rstdir=$(mktemp -d $TMPDIR/rstdir_XXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }   
        splittedfastapath=$rstdir/fasta
        $binpath/splitfasta.py -i $infile -nseq 1 -namemode 1 -outpath $splittedfastapath -ext fa
        mkdir -p $splittedfastapath
        numsplitted=`find $splittedfastapath -name "*.fa" | wc -l `
        for ((i=0;i<numsplitted;i+=1)); do
            echo "Running $i..."
            fafile=$splittedfastapath/${rootname}_${i}.fa
            if [ -s $fafile ] ; then 
                suboutdir=$rstdir/out$i
                mkdir -p $suboutdir
                RunPolyPhobius $fafile $outfile $suboutdir
                rm -rf $suboutdir
            fi 
        done
        /bin/rm -rf $rstdir
        echo "Finished at `/bin/date`"
        res2=$(/bin/date +%s.%N)
        printf "Running time for predicting %d sequences is: %.3F seconds\n" $numseq $(echo "$res2 - $res1"|/usr/bin/bc )
        echo "Result output to $outfile"
        echo "See log file at $logfile"
        if [ -s "$errfile" ]; then 
            echo "topcons exit with errors. Check the error log file '$errfile' for details"
        fi
    fi
fi

#clean TMPDIR
/bin/rm -rf $TMPDIR
#echo $TMPDIR
