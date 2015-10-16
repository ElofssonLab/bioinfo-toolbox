#!/bin/bash

# ChangeLog 2011-10-26 
#    topo2TMfrag.py and topoAddDGscore.py both are updated for that
#    1. accept topology with gaps
#    2. seqID must be supplied in the DG file, thus with three columns
#       in the format 
#       SeqID amino-acid-fragment dgValue

progname=`basename $0`
usage="
usage:  $progname TM2SEQSegmentFile -o [OUTFILE]

Given the TM2SEQ-segment file and append the DG values

  -h, --help     Print this help message and exit

Created 2013-04-18, updated 2013-04-18, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

AddAbsolutePath(){ #$path#{{{

    local var=$1
    if [ "${var:0:1}" != "/" ];then
        var=$PWD/$var # add the absolut path
    fi
    echo $var
    return 0
}
#}}}
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
TM2SEQSegmentAppendDGValue() { # $infile $outfile
    local infile=$1
    local outfile=$2
    tmpdir=$(mktemp -d /tmp/tmpdir.TM2SEQSegmentAppendDGvalue.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }
    trap 'rm -rf "$tmpdir"' INT TERM EXIT
    helixfile1=$tmpdir/helix1.txt
    helixfile2=$tmpdir/helix2.txt
    awk '/^[^#]/{print $(7)}' $infile > $helixfile1
    awk '/^[^#]/{print $(8)}' $infile > $helixfile2
    Nline=`grep -v "^#" $infile | wc -l`

    local dgcalcresult1=$tmpdir/dgcalcresult1.txt
    local dgcalcresult2=$tmpdir/dgcalcresult2.txt
    echo "$dgpredpath/calc_dG.pl $helixfile1  -outformat text > $dgcalcresult1"
    $dgpredpath/calc_dG.pl $helixfile1  -outformat text > $dgcalcresult1

    echo "$dgpredpath/calc_dG.pl $helixfile2  -outformat text > $dgcalcresult2"
    $dgpredpath/calc_dG.pl $helixfile2  -outformat text > $dgcalcresult2


    local dgscorelist1=$tmpdir/dgscorelist1.txt
    local dgscorelist2=$tmpdir/dgscorelist2.txt
    awk '{if ($1 != "Sequence" && $1 != "")print $2}' $dgcalcresult1 > $dgscorelist1
    awk '{if ($1 != "Sequence" && $1 != "")print $2}' $dgcalcresult2 > $dgscorelist2
    nline1=`cat $dgscorelist1 |wc -l`
    nline2=`cat $dgscorelist2 |wc -l`
    if [ "$nline1" != "$Nline" ]; then
        echo "Num1 ($nline1) != numRecord ($Nline)" >&2
    fi
    if [ "$nline2" != "$Nline" ]; then
        echo "Num2 ($nline2) != numRecord ($Nline)" >&2
    fi

    if [ "$outfile" != "" ]; then
        paste $infile $dgscorelist1 $dgscorelist2 > $outfile
    else
        paste $infile $dgscorelist1 $dgscorelist2
    fi

    rm -rf $tmpdir

}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

infile=
outfile=
binpath=`dirname $0`
dgpredpath=$binpath/../../../program/dgpred_standalone

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
            -o|--o) outfile=$2;shift;;
            -nc|--not-clean|-not-clean) isClean=0;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

if [ "$infile" == ""  ]; then
    echo "infile not set. exit" >&2
    exit
elif [ ! -s "$infile" ]; then
    echo "infile $infile does not exist or empty. exit" >&2
    exit
fi

IsProgExist  $dgpredpath/calc_dG.pl 

TM2SEQSegmentAppendDGValue $infile $outfile

