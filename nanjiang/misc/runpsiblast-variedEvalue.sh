#!/bin/bash
# 2009-06-11 
# run PSI-BLAST with varied E-value threshold for each item
# run psiblast for a given idListFile
# blastpgp version 2.2.17
# database=$DATADIR/blastdb/nr_20080117_pfilt
# iteration = 4
# threshold = varied

usage="
usage: ./runpsiblast-variedEvalue.sh [options] -l IDList  -pdbaa path
 no blastfile output
 Note            : IDList must have e-values
  -b int         : begin index for id, default = 0
  -e int         : end index for id, default = 0
  -a int         : num processor to run, default = 2
  -j int         : num of iterations for blastpgp, default=4
  --outpath path : output path for the blast result
  -db       file : set the blastdb, default=nr_20080117_pfilt
  -bb       int : __blastpgp_num_align_b
  -bm       int : __blastpgp_format_m
  -blastpgpopt  : blastpgp options
  -sleep <int>   : set sleep time in seconds, this is only for submitting batch jobs on multiple hosts
                 : to wait a certain amount of time, so that the login will not be blocked
                 : when the program starts
  --pdbaa|--aapath path   : the path to the sequence files in FASTA format

 begin is starting from 0, that is 
 b=0, e=1 will run on the first item 

 The low complexity filtering function is turned on
 Do not output the alignment file 
 The default evalue = 0.001

Updated 2010-08-27, Nanjiang Shu
"
function PrintHelp()
{
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

function RunPSIBLAST() # $id $e-value#{{{
#run psiblast , given the chain id and e-value threshold
{
    local id=$1
    local evalue=$2
    if  [ "$evalue" == "" ] ; then
        echo "Warning! line without E-value, set to 0.001"  >&2
        evalue=0.001
    fi

    local seqFile=$pdbaapath/$id.aa
    local chkfile=$outpath/$id.chk
    local pssmfile=$outpath/$id.pssm
    local blastfile=$outpath/$id.blast
    echo  "$BLASTBIN/blastpgp  -i $seqFile -h $evalue -d $database -j $iteration -a $numProcessor  -o $blastfile -Q $pssmfile -C $chkfile -b ${__blastpgp_num_align_b} -m ${__blastpgp_format_m} $blastpgpopt"
    $BLASTBIN/blastpgp -i $seqFile -h $evalue -d $database -j $iteration -a $numProcessor -o $blastfile -Q $pssmfile -C $chkfile -b ${__blastpgp_num_align_b} -m ${__blastpgp_format_m} $blastpgpopt
}
#}}}

BLAST_VERSION=blast-2.2.17
BLASTBIN=$HOME/usr/share/blast/$BLAST_VERSION/bin
BLASTMAT=$HOME/usr/share/blast/$BLAST_VERSION/data

nodename=`uname -n`
case $nodename in
    shu*) 
    BLASTDB=/data1/blastdb 
    BLASTBIN=$CASIODATA3/usr/share/blast/$BLAST_VERSION/bin
    BLASTMAT=$CASIODATA3/usr/share/blast/$BLAST_VERSION/data
    ;;
    casio*) 
    BLASTDB=/data/blastdb 
    BLASTBIN=$CASIODATA3/usr/share/blast/$BLAST_VERSION/bin
    BLASTMAT=$CASIODATA3/usr/share/blast/$BLAST_VERSION/data
    ;;
    scooby*) 
    BLASTDB=/do3/nanjiang/data/blastdb 
    BLASTBIN=$CASIODATA3/usr/share/blast/$BLAST_VERSION/bin
    BLASTMAT=$CASIODATA3/usr/share/blast/$BLAST_VERSION/data
    ;;
    hexa*|yarmulke*|pretoria*|linux*|grisman*)  #servers at cs.duke.edu
#    BLASTDB=/var/tmp/local/nanjiang/data/blastdb
    BLASTBIN=$HOME/usr/share/blast/$BLAST_VERSION/bin
    BLASTMAT=$HOME/usr/share/blast/$BLAST_VERSION/data
    ;;
    illergard*)
    BLASTDB=/nanjiang/data/blastdb
    BLASTBIN=$DATADIR1/usr/share/blast/$BLAST_VERSION/bin
    BLASTMAT=$DATADIR1/usr/share/blast/$BLAST_VERSION/data
    ;;
    *pdc.kth.se)
    BLASTDB=$DATADIR/blastdb
    BLASTBIN=$DATADIR1/usr/share/blast/$BLAST_VERSION/bin
    BLASTMAT=$DATADIR1/usr/share/blast/$BLAST_VERSION/data
    ;;
esac

#export BLASTDB
export BLASTBIN
export BLASTMAT
export BLAST_VERSION


begin=0
end=0

idtype=1
threshold=0.001
iteration=4
database=nr_20080117_pfilt
outpath=./variedEvalue
numProcessor=2
__blastpgp_num_align_b=100
__blastpgp_format_m=0
blastpgpopt=

sleepTime=0

# ##########to be modified each time
pdbaapath=
IDLIST=
# note: testVaridEvalue.list has two columns, one is the chain ID and the other is the E-value threshold
# #########to be modified each time

while [ "$1" != "" ]; do
	case $1 in
		-h|--help)PrintHelp;exit;;
		-b) begin=$2;shift ;;
		-e) end=$2;shift ;;
		-bb)__blastpgp_num_align_b=$2;shift ;;
		-bm)__blastpgp_format_m=$2;shift ;;
		-a) numProcessor=$2;shift ;;
		-j) iteration=$2;shift ;;
		--outpath|-outpath) outpath=$2;shift ;;
		--pdbaa|-pdbaa|--aapath|-aapath) pdbaapath=$2;shift ;;
		--blastpgpopt|-blastpgpopt) blastpgpopt="$2";shift ;;
		--db|-db) database=$2;shift ;;
        -sleep|--sleep) sleepTime=$2;shift;;
        -l) IDLIST=$2;shift;;
        *)  echo "Wrong argument $1"; exit;;
	esac
    # Shift all the parameters down by one
	shift
done

if [ ! -d $outpath ]; then
    mkdir -p $outpath
fi

if [ $sleepTime -ge 1 ]; then
    sleep ${sleepTime}s
fi

if [ ! -d $pdbaapath ]; then
    echo "Error! pdbaapath \"$pdbaapath\" does not exist!"
    exit
fi

if [ ! -f $IDLIST ]; then
    echo "Error! idListFile \"$IDLIST\" does not exist!"
    exit
fi

# read $FILE using the file descriptors
((cnt=0))
exec 3<&0
exec 0<$IDLIST
while read line
do
    # use $line variable to process line in processLine() function
    if [ $cnt -ge $begin -a $cnt -lt $end ]; then 
        echo   "$cnt  $line"
        RunPSIBLAST $line
        echo
    fi
    ((cnt++))
done
exec 0<&3

