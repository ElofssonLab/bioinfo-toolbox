#!/bin/bash

# build HHprofile 

# hhblits -i seqfile -oa3m seqfile.a3m -n 2 -d dbname
# hhmake -i seqfile.a3m -o outdir/seqfile.hmm
# ChangeLog 2012-07-06
#   option -no-overwrite added
# ChangeLog 2014-09-04
#   since hhsuite database after 2013 the untar file are located in a separate
#   folder, the location of the database name and *.tar.gz file should be named
#   accordingly
# ChangeLog 2014-09-08
#   1. the result file should be first write to a tmp folder and only when all
#   iterations are finished, then copied to the outpath
#   2. add the option -para, number of parallel hhblits jobs
#   3. add the option -cpu, number of cpus to use for each hhblits job
#   
#
nodename=`uname -n`
TMPDIR=/scratch/$$
case $nodename in 
    *pdc.kth.se) 
    datadir3=$DATADIR3 
    TMPDIR=/scratch/$$
    ;;
    tintin*|triolith*|*uppmax.uu.se)
    datadir3=$DATADIR3
    if [ $SNIC_TMP ]; then
        TMPDIR=$SNIC_TMP
    else
        if [ -d /scratch ]; then 
            TMPDIR=/scratch/$$
        else
            TMPDIR=/tmp/$$
        fi
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
if [ ! -d $TMPDIR ] ; then 
    mkdir -p $TMPDIR
fi

binpath=$datadir3/bin

usage="
USAGE:  buildHHprofile.sh fasta-seq-file [-outpath DIR]

Description: Build profiles for HHsearch, a3m and hhm files
  -outpath DIR      Set output path, (default: ./)
  -hhdb    STR      Basename of the HH database, 
                    (default: \$HHDB/uniprot20_02Sep11)
  -cpu     INT      Set the number of cpu to use for each hhblits job, (default: 2)
  -para    INT      Set the number of paralell jobs for hhblits, (defalt: 1)
  -no-overwrite     Do not overwrite the already existing file
  -hhblitsopt \"STR\" Set hhblits opt
                    e.g. -cpu INT -n 
  -q                Quiet mode
  -h|--help         Print this help message and exit

Created 2012-05-16, updated 2014-09-08, Nanjiang Shu 
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
exec_cmd(){
    echo "$1"
    eval "$1"
}
BuildHHprofile() { # $outdir#{{{
    local seqfile=$1
    if [ ! -s "$seqfile" ] ; then 
        echo "Error! seqfile '$seqfile' does not exit. Ignore." >&2
        return 1
    else
        local btname=`basename $seqfile`
        local rtname=${btname%.*}
        local tmpa3mfile=$seqfile.a3m
        local tmphhmfile=$seqfile.hhm

        local a3mfile=$outpath/$rtname.a3m
        local hhmfile=$outpath/$rtname.hhm

        if [ -s $a3mfile -a -s $hhmfile -a $isOverWrite -eq 0 ] ; then 
            echo "$a3mfile and $hhmfile exist. ignore"
        else
            exec_cmd "hhblits -i $seqfile -oa3m $tmpa3mfile -n 2 -d $hhdb -cpu $CPUperJob $hhblitsopt"
            exec_cmd "hhmake -i $tmpa3mfile -o $tmphhmfile"
            exec_cmd "/bin/cp -f $tmpa3mfile $a3mfile"
            exec_cmd "/bin/cp -f $tmphhmfile $hhmfile"
        fi
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
isOverWrite=1
outpath=./
infile=
hhdb=$HHDB/uniprot20_02Sep11
hhblitsopt=
CPUperJob=2
parallel=1

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
            -hhdb|--hhdb) hhdb=$2;shift;;
            -cpu|--cpu) CPUperJob=$2;shift;;
            -para|--para) parallel=$2;shift;;
            -hhblitsopt|--hhblitsopt) hhblitsopt="$2";shift;;
            -no-overwrite|--no-overwrite) isOverWrite=0;;
            -q|-quiet|--quiet) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

IsProgExist $binpath/splitfasta.py
IsProgExist $binpath/countseq.py

if [ ! -d $outpath ]; then 
    mkdir -p $outpath
fi
mkdir -p $outpath/log

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

    logfile=$outpath/log/${rootname}_buildHHprofile.log
    logfile=`readlink -f $logfile`

    if [ -s $logfile ]; then 
        cat  /dev/null > $logfile
    fi
    numseq=`$binpath/countseq.py $infile -nf `
    if [ "$numseq" == "" ]; then 
        echo "Fatal. countseq.py error. Exit."  >&2
        exit 1
    elif [ $numseq -le 0 ] ; then
        echo "Zero sequence found in file '$infile'. Exit."  >&2
        exit 1
    else
        tmpseqpath=$TMPDIR/aaseq
        echo "$numseq sequences are going to be processed by buildHHprofile"
        runningtime=`expr $numseq \* 2 `
        echo "It takes about $runningtime minutes to run."
        echo python $binpath/splitfasta.py -i $infile -ext aa -outpath $tmpseqpath
        python $binpath/splitfasta.py -i $infile -ext aa -outpath $tmpseqpath
        case $nodename in 
            tintin*|triolith*|*uppmax*)
                basename_hhdb=`basename $hhdb`
                dirname_hhdb=`dirname $hhdb`
                nametype=0
                mkdir -p $TMPDIR/hhsuite
                tarball=$hhdb.tar.gz
                if [ ! -f $tarball ];then
                    tarball=$dirname_hhdb.tar.gz #note for uniprot20_2013*, the tarball is at the folder one level up
                    nametype=1
                fi

                if [ ! -f $tarball ];then
                    echo "Database $tarball does not exist. exit $0"
                    exit
                fi
                cp $tarball $TMPDIR/hhsuite
                currdir=$PWD
                cd $TMPDIR/hhsuite
                tar -xvzf $tarball
                cd $currdir
                if [ "$nametype" == "0" ];then
                    hhdb=$TMPDIR/hhsuite/$basename_hhdb
                else
                    hhdb=$TMPDIR/hhsuite/$basename_hhdb/$basename_hhdb
                fi
                ;;
        esac

        ((i=0))
        for seqfile in $(find $tmpseqpath -name "*.aa"); do
            if [ ! -s "$seqfile" ]; then 
                echo seqfile $seqfile does not exist or empty, ignore.
                continue
            fi
            echo "Running $i..."
            ((i++))
            #echo "BuildHHprofile $seqfile &"
            BuildHHprofile $seqfile &
            if [[ $((i%$parallel)) -eq 0 ]];then
                wait
                echo count=$i, parallel=$parallel, CPUperJob=$CPUperJob, wait...
            fi
        done
        wait
        rm -rf $tmpseqpath
        echo "Finished at `/bin/date`"
        res2=$(/bin/date +%s.%N)
        printf "Running time for %d sequences is: %.3F seconds\n" $numseq $(echo "$res2 - $res1"|/usr/bin/bc )
        echo "See log file at $logfile"
    fi
fi

rm -rf $TMPDIR
