#!/bin/bash

usage="
usage  seqpairaln_to_topopaircmp.sh seq-pairaln-file -outpath DIR

Description: from seq-pairaln-file to pairwise topology comparison 

    -outpath DIR   Set dir for output result, (default: ./result)
    -topofile FILE Set topology file
    -seqdb    STR  Set sequence database file
    -dupfile FILE 
    -signalp FILE 

Created 2012-03-22, updated 2012-07-06, Nanjiang Shu 
"
SeqPairAln_To_TopoPairCmp() { #{{{
    local seqAlnFile=$1
    basename=`basename $seqAlnFile`
    rootname=${basename%.*}
    topowithdgscoreFile=$outpath/$rootname.topoaln.topowithdgscore

    if [ ! -s $topowithdgscoreFile ]; then 
        echo "Step 1. Get topology alignment from amino acid sequence alignment"
        topoAlnFile=$outpath/$rootname.topoaln.fa
        echo "$binpath/matchMSAtopo -msa $seqAlnFile -topo $topofile -o $topoAlnFile"
        $binpath/matchMSAtopo -msa $seqAlnFile -topo $topofile -o $topoAlnFile -ignore-badseq no

        echo "Step 2. Add DG values of all TM regions to the topoAlnFile"
        echo "$binpath/getDGvalueTMOfTopo.sh -topo $topoAlnFile -aa $dbfile -outpath $outpath"
        $binpath/getDGvalueTMOfTopo.sh -topo $topoAlnFile -aa $dbfile -outpath $outpath
    fi

    echo "Step 3. Compare pairwise topology alignments"
    cmpresultfile=$outpath/$rootname.paircmp
    echo "$binpath/compareMSATopo.py -mode 0 $topowithdgscoreFile -o $cmpresultfile -dupfile $dupfile -signalp $signalpfile"
    python $binpath/compareMSATopo.py -mode 0 $topowithdgscoreFile -o $cmpresultfile -dupfile $dupfile -signalp $signalpfile

    echo "Result output to $cmpresultfile"
}
#}}}

commandline="$0 $*"
rundir=`dirname $0`

binpath=$rundir
seqAlnFile=
dbfile=$rundir/../pfamAna/pairwise/all/pfamfullseq.selTM_uniq0.db0.db  
topofile=$rundir/../pfamAna/pairwise/all/predTM/pfamfullseq.selTM_uniq0.db0.db.topo
outpath=result
dupfile=
signalpfile=

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        seqAlnFile=$1
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) echo "$usage"; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -topofile|--topofile) topofile=$2;shift;;
            -dupfile|--dupfile) dupfile=$2;shift;;
            -signalp|--signalp) signalpfile=$2;shift;;
            -seqdb|--seqdb) dbfile=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -nc|-not-clean|--not-clean) isClean=0;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        seqAlnFile=$1
    fi
    shift
done

if [ "$seqAlnFile" == "" ]; then
    echo "seqAlnFile not set. Exit." >&2
    exit 1;
elif [ ! -s "$seqAlnFile" ] ; then 
    echo "seqAlnFile '$seqAlnFile' does not exist or empty. Exit." >&2
    exit 1;
fi

if [ ! -d "$outpath" ]; then 
    mkdir -p $outpath
fi

echo "Command line: $commandline"
echo "Start at    `date`"


SeqPairAln_To_TopoPairCmp $seqAlnFile

echo "Finished at `date`"

