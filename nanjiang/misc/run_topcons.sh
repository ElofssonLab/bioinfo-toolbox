#!/bin/bash

maxseq=1000

# Note 2012-04-10:
# topcons may not run successfully when the input sequence is more than
# 10

step=1

path_topcons=/server/var/www/topcons.cbr.su.se/docroot/topcons/
binpath=$DATADIR3/wk/MPTopo/src/

nodename=`uname -n`
TMPDIR=/scratch/$$
case $nodename in 
    *pdc.kth.se) 
    datadir3=$DATADIR3 
    TMPDIR=/scratch/$$
    path_topcons=$datadir3/share/var/www/topcons.cbr.su.se/docroot/topcons/
    ;;
    *uppmax.uu.se)
    datadir3=$DATADIR3
    if [ $SNIC_TMP ]; then
        TMPDIR=$SNIC_TMP
    else
        TMPDIR=/tmp/$$
    fi
    path_topcons=$datadir3/share/var/www/topcons.cbr.su.se/docroot/topcons/
    ;;
    trio*|*nsc.liu.se)
    datadir3=$DATADIR3
    if [ $SNIC_TMP ]; then
        TMPDIR=$SNIC_TMP
    else
        TMPDIR=/tmp/$$
    fi
    path_topcons=$datadir3/share/var/www/topcons.cbr.su.se/docroot/topcons/
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

if [ "$SNIC_SITE" == "nsc" ];then
    datadir3=$DATADIR3
    TMPDIR=$SNIC_TMP
fi

if [ ! -d $TMPDIR ] ; then 
    mkdir -p $TMPDIR
fi

#ChangeLog 2012-04-03 
#    methodstring added, before, I forgot to add replace the method string 
#ChangeLog 2012-06-07 
#    when running on uppmax, reset BLASTDB to TMPDIR/blastdb
#    and the script topcons/fa2prfs.sh is changed a bit
#ChangeLog 2014-10-30
#    option -topcons implemented
#ChangeLog 2015-01-08 
#   User set TMPDIR by -tmpdir


usage="
usage: run_topcons.sh fasta-seq-file [-outpath DIR]

Description: 
    Run web-server version of topcons locally
    output file will be $outpath/$rootname.topcons.allinfo

  -outpath DIR      Set output path, (default: ./)
  -topcons DIR      Set the path for topcons, 
                    (default: $path_topcons)
  -binpath DIR      Set path for scripts, (default: $binpath)
  -copyplot         Copy plot to outpath
  -tmpdir  DIR      User set TMPDIR
  -nc               Keep the temporary data
  -q                Quiet mode
  -h,--help         Print this help message and exit

Created 2012-04-28, updated 2015-01-08, Nanjiang Shu
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
RunTopcons(){ # $outdir#{{{
    local outdir=$1
    if [ ! -s "$outdir/query.tmp.fa" ] ; then 
        echo "Error! seqfile '$outdir/query' does not exit. Ignore." >&2
        return 1
    else
        outdir=`readlink -f $outdir`
        currdir=$PWD
        cd $path_topcons
        #./check_sequence.pl $outdir $maxseq 1> /dev/null 2>> $logfile
#        ./run_all.pl $outdir 1>$outdir/query.result.txt /dev/null 2>> $logfile
        dummyrstdir="dummyrstdir"
        dummyurl="dummyurl"

        echo "./run_all_html.pl $outdir $dummyrstdir $dummyurl $path_topcons 1> /dev/null 2>> $errfile"
        ./run_all_html.pl $outdir $dummyrstdir  $dummyurl $path_topcons 1> /dev/null 2>> $errfile
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
isClean=1
outpath=./
infile=
isCopyplot=0

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
            -topcons|--topcons) path_topcons=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -copyplot|--copyplot) isCopyplot=1;;
            -binpath|--binpath) binpath=$2;shift;;
            -tmpdir|--tmpdir) TMPDIR=$2/tmp_$$;shift;;
            -nc|--nc) isClean=0;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

if [ ! -d "$TMPDIR" ];then
    mkdir -p "$TMPDIR" ||  { echo "Failed to create temp dir \"$TMPDIR\"" >&2; exit 1; }
fi

echo TMPDIR=$TMPDIR
echo DATADIR3=$DATADIR3

IsProgExist $binpath/countseq.py
IsProgExist $binpath/catfasta.py
IsPathExist $path_topcons

/bin/mkdir -p $outpath

echo "path_topcons=$path_topcons"

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

    outfile=$outpath/$rootname.topcons.result.txt
    errfile=$outpath/$rootname.topcons.err
    logfile=$outpath/$rootname.topcons.log

    outfile=`readlink -f $outfile`
    errfile=`readlink -f $errfile`
    logfile=`readlink -f $logfile`

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
        echo "$numseq sequences are going to be predicted by topcons."
        runningtime=`expr $numseq / 2 \* 60`
        echo "It takes about $runningtime seconds to run the prediction."
        case $nodename in 
            *.pdc.kth.se|*uppmax.uu.se)
                tmpinfile=$(mktemp $TMPDIR/seqfile_XXXXXXX.fa) || { echo "Failed to create temp file" >&2; exit 1; }   
                /bin/cat $infile > $tmpinfile
                infile=$tmpinfile
                export BLASTDB=$TMPDIR/blastdb
                mkdir -p $BLASTDB
                echo "/bin/cp -f $path_topcons/uniref90.mem.fasta* $BLASTDB/"
                /bin/cp -f $path_topcons/uniref90.mem.fasta* $BLASTDB/
                ;;
        esac
        splittedpath=$TMPDIR/splittedseq
        mkdir -p $splittedpath
        $binpath/splitfasta.py -i $infile -nseq 1 -outpath $splittedpath -ext fa -nameseq
        numsplitted=`find $splittedpath -name "*.fa" | wc -l `
        for ((i=0;i<numsplitted;i+=step)); do
            splittedfastafile=$splittedpath/${rootname}_${i}.fa
            rstdir=$(mktemp -d $TMPDIR/rstdir_XXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }   
            cat $splittedfastafile > $rstdir/query.tmp.fa
            echo "Running $i..."
            RunTopcons $rstdir
            if [ -s "$rstdir/query.txt" ] ; then
                cat $rstdir/query.txt >> $outfile
            fi
            if [ $isCopyplot -eq 1 ];then
                /bin/cp -rf $rstdir $outpath/${rootname}_${i}
                out_result_path=$outpath/${rootname}_${i}/result/dummyrstdir/
                mkdir -p $out_result_path
                declare -a sourcefilelist=("query.png"  "query.large.png"  "query.topcons.png"  "query.topcons.large.png"  "query.txt"  "query.blast")
                declare -a targetfilelist=("result.png"  "result.large.png"  "result.topcons.png"  "result.topcons.large.png"  "topcons.txt"  "blast.txt")
                numfile_to_copy=${#sourcefilelist[@]}
                for ((j=0;j<numfile_to_copy;j++)); do
                    /bin/cp -f $rstdir/${sourcefilelist[$j]} $out_result_path/${targetfilelist[$j]}
                done

            else
                if [ $isClean -eq 1 ];then
                    /bin/rm -rf $rstdir
                fi
            fi
        done
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
if [ $isClean -eq 1 ];then
    /bin/rm -rf $TMPDIR
else
    echo "Temporary data are kept at $TMPDIR"
fi
#echo $TMPDIR
