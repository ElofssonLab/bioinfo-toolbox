#!/bin/bash


progname=`basename $0`
usage="
usage:  $progname pairinfofile -seqdb dbname -topodb dbname -o [OUTFILE]

Description:
    Given pairinfo and append krbias

  -h, --help     Print this help message and exit

Created 2013-07-01, updated 2013-07-01, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

exec_cmd(){
    echo "$1"
    eval "$1"
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
PairInfoAppendKRbias() { # $infile $outfile#{{{
    local infile=$1
    local outfile=$2
    local tmpdir=$(mktemp -d /tmp/tmpdir.pairinfoAppendKRbias.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }
    trap 'rm -rf "$tmpdir"' INT TERM EXIT
    local idlistfile=$tmpdir/idlistfile.txt
    local topofile=$tmpdir/topofile.txt
    local seqfile=$tmpdir/seqfile.txt
    awk '/^[^#]/{print $1; print $2;}' $infile > $idlistfile
    exec_cmd "$DATADIR3/wk/MPTopo/src/my_extractdb.py -dbname $topodb -l $idlistfile -o $topofile"
    exec_cmd "$DATADIR3/wk/MPTopo/src/my_extractdb.py -dbname $seqdb -l $idlistfile -o $seqfile"

    exec_cmd "$DATADIR3/wk/MPTopo/src/test_calculate_KR_bias.py -maxdist 12 -flankwin 0 -topofile $topofile -seqfile $seqfile -outpath $tmpdir"

    local krlistfile=$tmpdir/seqfile.krlist.txt

    nline1=`cat $idlistfile |wc -l`
    nline2=`cat $krlistfile |wc -l`
    if [ "$nline1" != "$nline2" ]; then
        echo "numID ($nline1) != numKRbias ($nline2)" >&2
    fi

    local krbiasfile=$tmpdir/seqfile.krbias
    awk '{if(NR%2==1){printf("%s",$3)}else{printf("\t%s\n",$3)}}' $krlistfile | sed 's/[():]//g'  > $krbiasfile
    if [ "$outfile" != "" ]; then
        paste $infile $krbiasfile > $outfile
    else
        paste $infile $krbiasfile
    fi

    exec_cmd "/bin/cp -f $krlistfile ."


}
#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

infile=
outfile=
binpath=`dirname $0`
topodb=$DATADIR3/wk/MPTopo/pfamAna_refpro/pred_topcons/refpro20120604-celluar.selmaxlength-m1.topcons.result_TOPCONS.m1.agree-55
seqdb=$DATADIR3/wk/MPTopo/pfamAna_refpro/cellular_filter_fragment/pairwise/withinPfam/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20.nr100.filter.fragmented.uniq.pfam
isTopoDBSet=0

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
        -seqdb|--seqdb) seqdb=$2;shift;;
        -topodb|--topodb) topodb=$2; isTopoDBSet=1;shift;;
            -nc|--not-clean|-not-clean) isClean=0;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

if [ $isTopoDBSet -eq 0 ];then
    case "$PWD/$infile" in 
        *topcons_single*)
            topodb=$DATADIR3/wk/MPTopo/pfamAna_refpro/pred_topcons_single_method4/refpro20120604-celluar.selmaxlength-m1.topcons-single_topcons_single.m1.agree-44
            ;;
        *topcons*)
            topodb=$DATADIR3/wk/MPTopo/pfamAna_refpro/pred_topcons/refpro20120604-celluar.selmaxlength-m1.topcons.result_TOPCONS.m1.agree-55
            ;;
    esac
fi

if [ "$infile" == ""  ]; then
    echo "infile not set. exit" >&2
    exit
elif [ ! -s "$infile" ]; then
    echo "infile $infile does not exist or empty. exit" >&2
    exit
fi


PairInfoAppendKRbias $infile $outfile

