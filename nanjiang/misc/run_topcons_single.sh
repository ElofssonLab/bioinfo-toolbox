#!/bin/bash

maxseq=1000
# Note 2012-04-10:
# topcons_single may not run successfully when the input sequence is more than
# 10
step=5
path_topcons_single=$DATADIR3/share/var/www/single.topcons.net/docroot/topcons-single/
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
    datadir3=$DATADIR3
    if [ -d /scratch ]; then 
        TMPDIR=/scratch/$$
    else
        TMPDIR=/tmp/$$
    fi
    ;;
esac

#ChangeLog 2012-04-03 
#    methodstring added, before, I forgot to add replace the method string 
#ChangeLog 2012-05-24
#   options added
#    -beg INT
#    -end INT 

export HMMTOP_ARCH=$datadir3/share/hmmtop/bin/hmmtop.arch
export HMMTOP_PSV=$datadir3/share/hmmtop/bin/hmmtop.psv
export MEMSAT_DIR=$datadir3/share/memsat/
path_topcons_single=$datadir3/share/var/www/single.topcons.net/docroot/topcons-single/
binpath=$datadir3/bin

usage="
usage:   run_topcons_single.sh fasta-seq-file [-outpath DIR]
Description: Run web-server version of topcons-single locally
             output file will be $outpath/$rootname.topcons-single.allinfo
  -outpath DIR      Set output path, (default: ./)
  -continue|-cont   Write after the previous result file
  -beg     INT      Index of the first sequence to be run, (default: 0)
  -end     INT      Index of the last sequence to be run, (default: numseq)
  -method  INT      Set the number of methods to use, (default: 4)
                    3: SCAMPI_single, HMMTOP, MEMSAT-1.0
                    4: SCAMPI_single, HMMTOP, S-TMHMM, MEMSAT-1.0
                    5: SCAMPI_single, HMMTOP, TopPred, S-TMHMM, MEMSAT-1.0
                    6: SCAMPI_single, HMMTOP, TopPred, Phobius, S-TMHMM, MEMSAT-1.0
  -q                Quiet mode
  -h|--help         Print this help message and exit
Created 2011-11-12, updated 2012-05-24, Nanjiang Shu 
"
function PrintHelp() {
    echo "$usage"
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
function RunTopconsSingle() { # $outdir#{{{
    local outdir=$1
    if [ ! -s "$outdir/query.tmp.fa" ] ; then 
        echo "Error! seqfile '$outdir/query.tmp.fa' does not exit. Ignore." >&2
        return 1
    else
        outdir=`readlink -f $outdir`
        currdir=$PWD
        cd $path_topcons_single
        ./check_sequence.pl $outdir $maxseq 1> /dev/null 2>> $errfile
        ./run_all.pl $outdir/ $methodstring 1> /dev/null 2>> $errfile
        cd $currdir
        return 0
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
method=4;
methodstring=stmhmm_hmmtop_memsat
g_start=0
g_end=2147483647
isContinue=0

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
            -beg|--beg) g_start=$2;shift;;
            -end|--end) g_end=$2;shift;;
            -method|--method) method=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -cont|--cont|--continue|-continue) isContinue=1;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

IsProgExist $binpath/countseq.py
IsProgExist $binpath/catfasta.py
IsPathExist $path_topcons_single

case $method in 
    3) methodstring=hmmtop_memsat ;;
    4) methodstring=hmmtop_stmhmm_memsat ;;
    5) methodstring=hmmtop_toppred_stmhmm_memsat ;;
    6) methodstring=hmmtop_toppred_phobius-seq_stmhmm_memsat ;;
    *) echo "Wrong method '$method'. Exit.";exit 1;;
esac


if [ "$infile" == ""  ]; then 
    echo "infile not set. Exit" >&2
    exit 1
elif [ ! -s "$infile" ] ; then
    echo "infile '$infile' not exist or empty. Exit. " >&2
else 
    if [ ! -d $TMPDIR ] ; then 
        mkdir -p $TMPDIR
    fi
    if [ ! -d $outpath ]; then
        mkdir -p $outpath
    fi
    echo "Command line: $commandline"
    echo "Start at    `/bin/date`"
    res1=$(/bin/date +%s.%N)
    basename=`basename $infile`
    rootname=${basename%.*}
    outfile=$outpath/$rootname.topcons-single.allinfo
    errfile=$outpath/$rootname.topcons-single.err

    outfile=`readlink -f $outfile`
    errfile=`readlink -f $errfile`

    if [ $isContinue -eq 0 ] ; then 
        cat  /dev/null > $outfile
        cat  /dev/null > $errfile
    fi
    numseq=`$binpath/countseq.py $infile -nf `
    if [ "$numseq" == "" ]; then 
        echo "Fatal. countseq.py error. Exit."  >&2
        rm -rf $TMPDIR
        exit 1
    elif [ $numseq -le 0 ] ; then
        echo "Zero sequence found in file '$infile'. Exit."  >&2
        rm -rf $TMPDIR
        exit 1
    else
        if [ $numseq -lt $g_end ]; then
            g_end=$numseq
        fi
        echo "$numseq sequences are going to be predicted by topcons-single."
        runningtime=`expr $numseq / 2 `
        echo "It takes about $runningtime seconds to run the prediction."
        case $nodename in 
            *.pdc.kth.se|*uppmax.uu.se)
                tmpinfile=$(mktemp $TMPDIR/seqfile_XXXXXXX.fa) || { echo "Failed to create temp file" >&2; rm -rf $TMPDIR; exit 1; }   
                /bin/cat $infile > $tmpinfile
                infile=$tmpinfile
                ;;
        esac
        for ((i=$g_start ; i<=$g_end;i+=step)); do
            rstdir=$(mktemp -d $TMPDIR/rstdir_XXXXXXX) || { echo "Failed to create temp dir" >&2; rm -rf $TMPDIR; exit 1; }   
            ((beg=i))
            ((end=i+step))
            $binpath/catfasta.py -b $beg -e $end $infile > $rstdir/query.tmp.fa
            if [ -s "$rstdir/query.tmp.fa" ]; then 
                echo "Running $i..."
                RunTopconsSingle $rstdir
                if [ -s "$rstdir/all_info.txt" ] ; then
                    cat $rstdir/all_info.txt >> $outfile
                fi
            fi
            rm -rf $rstdir
        done
        echo "Finished at `/bin/date`"
        res2=$(/bin/date +%s.%N)
        printf "Running time for predicting %d sequences is: %.3F seconds\n" $numseq $(echo "$res2 - $res1"|/usr/bin/bc )
        echo "Result output to $outfile"
        if [ -s "$errfile" ]; then 
            echo "topcons-single exit with errors. Check the error log file '$errfile' for details"
        fi
    fi
fi

#clean TMPDIR
rm -rf $TMPDIR
