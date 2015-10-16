#!/bin/bash

binpath=/misc/illergarddata/wk/MPTopo/bin/
dgpredpath=/misc/illergarddata/program/dgpred/dgpred_standalone/
usage="
Usage:  selectDIFFTopo.sh seqID -binpath binpath -dgpredpath dgpredpath -datapath datapath

select DIFF topo homologues predicted with high probability

Options:
    -binpath    <dir> : set the path for executables, default=$binpath 
    -dgpredpath <dir> : set the path for executables or dgpred, default=$dgpredpath
    -datapath   <dir> : set the path storing the data
    -h|--help         : print this help message and exit

Created 2010-09-07, updated 2010-09-07, Nanjiang
nanjiang.shu@gmail.com
"
function PrintHelp()
{
    echo "$usage"
}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi


id=

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -binpath|--binpath) binpath=$2;shift;;
            -dgpredpath|--dgpredpath) dgpredpath=$2;shift;;
            -datapath|--datapath) datapath=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        id=$1
    fi
    shift
done

if [ "$id"  == "" ]; then 
    echo "Error! id not set. Exit..."
    exit
fi
if [ ! -d "$binpath" ]; then 
    echo "Error!  binpath= \"$binpath\" does not exist. Exit..."
    exit
fi
if [ ! -d "$dgpredpath" ]; then 
    echo "Error!  dgpredpath= \"$dgpredpath\" does not exist. Exit..."
    exit
fi
if [ ! -d "$datapath" ]; then 
    echo "Error!  datapath= \"$datapath\" does not exist. Exit..."
    exit
fi
# extract the amino acid fragment for TM regions
$binpath/topo2TMfrag.py  $datapath/$id.blast.topo -f $datapath/$id.blast.tmp.fa -o $datapath/$id.blast.tmp.helixfrag
$binpath/topo2TMfrag.py  $datapath/$id.topo -f $datapath/$id.fa > $datapath/$id.helixfrag

# calculate DG values for all TM segments
$dgpredpath/calc_dG.pl $datapath/$id.blast.tmp.helixfrag | html2text | perl -i.bk -pe 's/[^[:ascii:]]//g;' > $datapath/$id.blast.tmp.dg
$dgpredpath/calc_dG.pl $datapath/$id.helixfrag | html2text | perl -i.bk -pe 's/[^[:ascii:]]//g;' > $datapath/$id.dg

# add DG values to topo files
$binpath/topoAddDGscore.py $datapath/$id.topo -dg $datapath/$id.dg --outpath $datapath
$binpath/topoAddDGscore.py $datapath/$id.blast.topo -dg $datapath/$id.blast.tmp.dg --outpath $datapath

# add DG values to topology comparison files
$binpath/resultAddAlnDGScore.py -topocmp $datapath/$id.blast.topocmp_and_hmmscore -topoaln $datapath/$id.topoaln  -qdg $datapath/$id.dgscorelist -tdg $datapath/$id.blast.dgscorelist -o $datapath/$id.topocmp_and_hmmscore.topoalnwithdgscore

# select topology alignments based on DG values and HMM scores and PIDs and Evalues
$binpath/selTopoAlnWithDGScore.py $datapath/$id.topocmp_and_hmmscore.topoalnwithdgscore -o $datapath/$id.sel.topoalnwithdgscore
