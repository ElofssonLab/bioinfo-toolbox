#!/bin/bash
# filename: batall_pairwise_hhalign1.sh
# this script calls the script  $DATADIR3/wk/MPTopo/src/bat_pfamana_pairwise_hhalign.sh 

PrintHelp(){ #{{{
    echo "$usage"
}
#}}}
exec_cmd() { #{{{
    echo "$*"
    eval "$*"
} #}}}

#ChangeLog 2014-05-15
#   -subdir STR    setting subdirs to work with, if not set, the default subdir
#   list will be used
subdirlist_default="
topcons_single_agreement_44
topcons_single_agreement_43
topcons_single_agreement_42
topcons_single_agreement_41
topcons_single_agreement_40
topcons_agreement_55
topcons_agreement_54
topcons_agreement_53
topcons_agreement_52
topcons_agreement_51
topcons_agreement_50
"


progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname [-workdir DIR]
Options:
  -coverage FLOAT   coverage threshold for hhalign alignment, (default: 0.0)
  -tax      STR     Taxonomy group. all, bac, arc, pro, euk, (default: all)
                      all: Everything
                      bac: bacteria
                      arc: Archaea
                      pro: Prokaryote
                      euk: Eukaryote
  -pfamtype STR     Type of Pfam family, all, fam, dom, (default: Family)
                      all: Everything
                      fam: Family
                      dom: Domain
  -topotype STR     topology type, (default: all)
                      all: Everything
                      delsigp: proteins with signal peptide are removed
                      rmsp: signal peptide is masked
  -seqtype STR     sequence type, (default: all)
                      all: Everything
                      sw: swissprot
  -mp INT          pairwise comparion method, 0 or 1, (default: 3)


  -stopbeforehhalign  Whether stop before running hhalign, just create 
                      pair list
  -subdir  STR      Add subdirs to proceed, e.g. topcons_single_agreement_44
  -h, --help        Print this help message and exit

Created 2013-08-27, updated 2014-05-15, Nanjiang Shu

Note: subdirlist_default is 
$subdirlist_default
"
maxpairlist="
100
"

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

workdir=
coverage=0
pfamtype=Family
taxgroup=all
seqtype=all
topotype=all
isStopBeforeHHalign=0
pairwise_comparison_method=3
subdirlist=

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        echo Error! Wrong argument: $1 >&2; exit;
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -workdir|--workdir) workdir=$2;shift;;
            -coverage|--coverage) coverage=$2;shift;;
            -mact|--mact) mact=$2;shift;;
            -tax|--tax) taxgroup=$2;shift;;
            -pfamtype|--pfamtype) pfamtype=$2;shift;;
            -topotype|--topotype) topotype=$2;shift;;
            -seqtype|--seqtype) seqtype=$2;shift;;
            -subdir|--subdir) subdirlist="$subdirlist $2";shift;;
            -stopbeforehhalign|--stopbeforehhalign)isStopBeforeHHalign=1;;
        -mp|--mp) pairwise_comparison_method=$2;shift;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        echo Error! Wrong argument: $1 >&2; exit;
    fi
    shift
done

if [ "$workdir" == "" ];then 
    echo "workdir not set. exit"
    exit 1
fi

taxgroup=`echo $taxgroup | tr '[A-Z]' '[a-z]'`
case $taxgroup in 
    al*) taxgroup=all;;
    ar*) taxgroup=Archaea;;
    ba*) taxgroup=Bacteria;;
    pr*) taxgroup=Prokaryota;;
    eu*) taxgroup=Eukaryota;;
    *) echo "Error! Wrong taxgroup '$taxgroup'"; exit 1;;
esac

pfamtype=`echo $pfamtype | tr '[A-Z]' '[a-z]'`
case $pfamtype in 
    al*) pfamtype=all;;
    fa*) pfamtype=Family;;
    do*) pfamtype=Domain;;
    *) echo "Error! Wrong pfamtype '$pfamtype'"; exit 1;;
esac

topotype=`echo $topotype | tr '[A-Z]' '[a-z]'`
topotype_addname=
case $topotype in 
    al*) topotype=all; topotype_addname=;;
    del*) topotype=delSPseq; topotype_addname=.delSPseq;;
    rm*) topotype=RMSP; topotype_addname=.RMSP;;
    *) echo "Error! Wrong topotype '$topotype'"; exit 1;;
esac

seqtype=`echo $seqtype | tr '[A-Z]' '[a-z]'`
seqtype_addname=
case $seqtype in 
    al*) seqtype=all; seqtype_addname=;;
    sw*) seqtype=swissprot; seqtype_addname=.swissprot;;
    *) echo "Error! Wrong seqtype '$seqtype'"; exit 1;;
esac


opt=
if [ $isStopBeforeHHalign -eq 1 ];then
    opt="-stopbeforehhalign"
fi


if [ ! -d $workdir ];then
    mkdir -p $workdir
fi
currdir=$PWD

cd $workdir
echo "Now in the workdir: $workdir"

if [ ! -d log ];then
    mkdir -p log
fi

maxjob=4
count=0

if [ "$subdirlist" == "" ];then
    subdirlist="$subdirlist_default"
fi

echo "subdirlist to proceed is"
for subdir in $subdirlist; do
    echo -e "\t$subdir"
done

for subdir in $subdirlist; do
    for maxpair in $maxpairlist;do
        logfile=log/batall_hhalign.${subdir}.maxpair${maxpair}.log
        case $workdir in 
            *pfamAna_refpro*)
                exec_cmd "bash $DATADIR3/wk/MPTopo/src/bat_pfamana_pairwise_hhalign.sh -workdir $subdir -coverage $coverage -tax $taxgroup -pfamtype $pfamtype -topotype $topotype -seqtype $seqtype -maxpair $maxpair -mp $pairwise_comparison_method $opt > $logfile 2>&1 &"
                ;;
            *pfamAna_swissprot*)
                exec_cmd "bash $DATADIR3/wk/MPTopo/src/bat_pfamana_swissprot_pairwise_hhalign.sh -workdir $subdir -coverage $coverage -tax $taxgroup -pfamtype $pfamtype -topotype $topotype -seqtype $seqtype -maxpair $maxpair -mp $pairwise_comparison_method $opt > $logfile 2>&1 &"
                ;;
        esac
        ((count++))
        if [[ $((count%maxjob)) -eq 0 ]] ; then
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $count wait $subdir $maxpair
        fi
    done
done

wait

cat `find . -name "*.maxpair$maxpair.pairlist" | grep "topcons"`  | sort -u > dumped.maxpair$maxpair.pairlist

cd $currdir
