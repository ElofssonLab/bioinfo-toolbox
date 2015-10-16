#!/bin/bash

# ChangeLog 2011-10-26 
#    topo2TMfrag.py and topoAddDGscore.py both are updated for that
#    1. accept topology with gaps
#    2. seqID must be supplied in the DG file, thus with three columns
#       in the format 
#       SeqID amino-acid-fragment dgValue

usage="
usage:  getDGvalueTMOfTopo.sh -topo topoFile -aa aaSeqFile

Get the DG value of all TM regions given the topology file and the
corresponding amino acid sequence file.
Topology file can be with gaps
Both should be in the Fasta format 

The output file will be $outpath/$rootname(topoFile).topowithdgscore
if topo are not included, then the output file will be named as $outpath/$rootname(topoFile).dgtm

  -outpath DIR   Set output path
  -l       FILE  Set the fileListFile, a pair per line
  -q             Quiet mode
  -verbose       Print verbose information
  -h|--help      Print this help message and exit

Created 2011-09-08, updated 2012-03-09, Nanjiang Shu
"
function PrintHelp() {
    echo "$usage"
}

function AddAbsolutePath() #$path#{{{
{
    local var=$1
    if [ "${var:0:1}" != "/" ];then
        var=$PWD/$var # add the absolut path
    fi
    echo $var
    return 0
}
#}}}
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
function GetDGvalueTMofTopo() #$topoFile $aaSeqFile#{{{
{
# the most time consuming part is calc_dG.pl

    local topoFile=$1
    local aaSeqFile=$2
    local tmpdir=$(mktemp -d /tmp/tmpdir.getDGvalueTMOfTopo.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }   
    local basename=`basename "$topoFile"`
    local rootname=${basename%.*}

    local helixfragwithid=$tmpdir/helixfragwithid.txt
    if [ $verbose -eq 1 ]; then 
        echo "python $binpath/topo2TMfrag.py  $topoFile -f $aaSeqFile -printid yes -o $helixfragwithid" 
    fi
#     res1=$(/bin/date +%s.%N)
    python $binpath/topo2TMfrag.py  $topoFile -f $aaSeqFile -printid yes -o $helixfragwithid 
#     res2=$(/bin/date +%s.%N)
#     printf "Running time for %s: %.3F\n" "topo2TMfrag" $(echo "$res2 - $res1"|/usr/bin/bc )

    local helixfragwithoutid=$tmpdir/helixfragwithoutid.txt
    awk '{print $2}' $helixfragwithid > $helixfragwithoutid

    local dgcalcresult=$tmpdir/dgcalcresult.txt
    if [ $verbose -eq 1 ]; then 
        echo "$dgpredpath/calc_dG.pl $helixfragwithoutid  -outformat text > $dgcalcresult"
    fi
#     res1=$(/bin/date +%s.%N)
    $dgpredpath/calc_dG.pl $helixfragwithoutid  -outformat text > $dgcalcresult
#     res2=$(/bin/date +%s.%N)
#     printf "Running time for %s: %.3F\n" "calc_dG" $(echo "$res2 - $res1"|/usr/bin/bc )

    local dgscorelist=$tmpdir/dgscorelist.txt
    awk '{if ($1 != "Sequence" && $1 != "")print $2}' $dgcalcresult > $dgscorelist

    local dgscorelist_with_id_seq=$outpath/${rootname}.dgscorelist
    paste $helixfragwithid $dgscorelist > $dgscorelist_with_id_seq
    echo "$dgscorelist_with_id_seq output"
      
    if [ $verbose -eq 1 ]; then 
        echo "python $binpath/topoAddDGscore.py $topoFile -dg $dgscorelist_with_id_seq --outpath $outpath"
    fi
#     res1=$(/bin/date +%s.%N)
    python $binpath/topoAddDGscore.py $topoFile -dg $dgscorelist_with_id_seq --outpath $outpath
#     res2=$(/bin/date +%s.%N)
#     printf "Running time for %s:  %.3F\n" "topoAddDGscore" $(echo "$res2 - $res1"|/usr/bin/bc )
    if [ $isClean -eq 1 ] ;then
        rm -rf $tmpdir
    else
        echo "tmpdir saved at $tmpdir"
    fi
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=0
verbose=0
isClean=1
outpath=./
fileListFile=
topoFile=
aaSeqFile=

binpath=`dirname $0`
dgpredpath=$binpath/../../../program/dgpred_standalone

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -topo|--topo) topoFile=$2;shift;;
            -aa|--aa) aaSeqFile=$2;shift;;
            -l|--l|-list|--list) fileListFile=$2;shift;;
            -q|-quiet|--quiet) isQuiet=1;;
            -verbose|--verbose) verbose=1;;
            -nc|--not-clean|-not-clean) isClean=0;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        echo "Error! Wrong argument: $1">&2; exit;
    fi
    shift
done

isFileNotSet=1;
if [ "$topoFile" != "" -a "$aaSeqFile" != "" ] ; then
    isFileNotSet=0;
fi
if [ "$fileListFile" != ""  ]; then
    isFileNotSet=0;
fi
if [ $isFileNotSet -eq 1 ]; then 
    echo "$0: Error, input not set!" >&2
    exit
fi
mkdir -p $outpath
IsProgExist $binpath/topoAddDGscore.py
IsProgExist  $dgpredpath/calc_dG.pl 
IsProgExist  $binpath/topo2TMfrag.py

if [ "$topoFile" != "" -a "$aaSeqFile" != "" ]; then  
    GetDGvalueTMofTopo $topoFile $aaSeqFile
fi

if [ "$fileListFile" != ""  ]; then
    if [ -s "$fileListFile" ] ; then
        exec 3<&0
        exec 0<$fileListFile
        while read line
        do
            GetDGvalueTMofTopo $line
        done
        exec 0<&3
    else 
        echo "list file $fileListFile does not exist" >&2
    fi
fi

