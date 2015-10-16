#!/bin/bash
exec_cmd(){
    echo "$*"
    eval "$*"
}
rundir=`dirname $0`

# Filename: rmsigp.sh
# Description: process prediction topology by TOPCONS or TOPCONS-single and
#              also remove signal peptide prediction
# Author: Nanjiang Shu (nanjiang.shu@scilifelab.se)

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname -path_predtopo DIR -path_signalp -basename STR [-mode tp|tps]
Options:
  -path_predtopo DIR   Path for predicted topology
  -path_signalp DIR    Path for predicted signal peptide by signalp4
  -basename     STR    basename of the file
                       TOPCONS prediction will be \$path_predtopo/\$basename.topcons.result.txt
                       TOPCONS-single prediction will be \$path_predtopo/\$basename.topcons-single.allinfo
  -h, --help        Print this help message and exit

Created 2014-10-09, updated 2014-10-09, Nanjiang Shu 
"
PrintHelp(){ #{{{
    echo "$usage"
}
#}}}
ProcessPredTopo_tps(){ #{{{
    local predfile=$path_predtopo/${basename}.topcons-single.allinfo
    if [ ! -s "$predfile" ];then
        echo "$predfile does not exist or empty, ignore"
        return 1
    fi

# extract 
    exec_cmd "$binpath/topcons_single_to_fasta.py $predfile -outpath $path_predtopo -m 1"
    exec_cmd "$binpath/topcons_single_to_fasta.py $predfile -outpath $path_predtopo -m 0"
# select agreement
    local num_predictor=4
    local agreementfile=$path_predtopo/${basename}.topcons-single.agreement.stat.txt
    local topodb=$path_predtopo/${basename}.topcons-single_topcons_single.topo
    local topofile=
    local outfile=
    local tmproidlist=
    if [ ! -f  ${topodb}0.db ];then
        exec_cmd "$binpath/indexfasta.py $topodb"
    fi

    local numIDTList=`seq 0 $num_predictor`
    for numIDT in $numIDTList; do
        tmproidlist=$path_predtopo/${basename}_topcons_single.m1.agree-${num_predictor}${numIDT}.tmproidlist

        awk -v NIDT=$numIDT -v NPD=$num_predictor '
        {if($3==NPD && $2==NIDT)print $1}
        ' $agreementfile > $tmproidlist

        # extract topology
        topofile=$path_predtopo/${basename}_topcons_single.m1.agree-${num_predictor}${numIDT}.topo
        rm -f $topofile
        exec_cmd "$binpath/my_extractdb.py -dbname $topodb  -l $tmproidlist -o $topofile"
        exec_cmd "$binpath/indexfasta.py $topofile"
    done

    # 3. filter signal peptide

    for item in all 40 41 42 43 44; do
        case $item in 
            all)
                topofile=$path_predtopo/${basename}.topcons-single_topcons_single.topo
                outfile=$path_predtopo/${basename}.topcons-single_topcons_single.RMSP.topo
                ;;
            *)
                topofile=$path_predtopo/${basename}_topcons_single.m1.agree-$item.topo
                outfile=$path_predtopo/${basename}_topcons_single.m1.agree-$item.RMSP.topo
                ;;
        esac
        rm -f $outfile
        exec_cmd "$binpath/filterSignalPeptide.py -topo  $topofile -sig  $signalpfile -o $outfile"
        exec_cmd "$binpath/indexfasta.py $outfile"
    done
}
#}}}
ProcessPredTopo_tp(){ #{{{
    local predfile=$path_predtopo/${basename}.topcons.result.txt
    if [ ! -s "$predfile" ];then
        echo "$predfile does not exist or empty, ignore"
        return 1
    fi

# extract 
    exec_cmd "$binpath/topcons_to_fasta.py $predfile -outpath $path_predtopo -m 1"
    exec_cmd "$binpath/topcons_to_fasta.py $predfile -outpath $path_predtopo -m 0"
# select agreement
    local num_predictor=5
    local agreementfile=$path_predtopo/${basename}.topcons.result.agreement.stat.txt
    local topodb=$path_predtopo/${basename}.topcons.result_TOPCONS.topo
    local topofile=
    local outfile=
    local tmproidlist=
    if [ ! -f  ${topodb}0.db ];then
        exec_cmd "$binpath/indexfasta.py $topodb"
    fi

    local numIDTList=`seq 0 $num_predictor`
    for numIDT in $numIDTList; do
        tmproidlist=$path_predtopo/${basename}.topcons.result_TOPCONS.m1.agree-${num_predictor}${numIDT}.tmproidlist

        awk -v NIDT=$numIDT -v NPD=$num_predictor '
        {if($3==NPD && $2==NIDT)print $1}
        ' $agreementfile > $tmproidlist

        # extract topology
        topofile=$path_predtopo/${basename}.topcons.result_TOPCONS.m1.agree-${num_predictor}${numIDT}.topo
        rm -f $topofile
        exec_cmd "$binpath/my_extractdb.py -dbname $topodb  -l $tmproidlist -o $topofile"
        exec_cmd "$binpath/indexfasta.py $topofile"
    done

    # 3. filter signal peptide

    for item in all 50 51 52 53 54 55; do
        case $item in 
            all)
                topofile=$path_predtopo/${basename}.topcons.result_TOPCONS.topo
                outfile=$path_predtopo/${basename}.topcons.result_TOPCONS.RMSP.topo
                ;;
            *)
                topofile=$path_predtopo/${basename}.topcons.result_TOPCONS.m1.agree-$item.topo
                outfile=$path_predtopo/${basename}.topcons.result_TOPCONS.m1.agree-$item.RMSP.topo
                ;;
        esac
        rm -f $outfile
        exec_cmd "$binpath/filterSignalPeptide.py -topo  $topofile -sig  $signalpfile -o $outfile"
        exec_cmd "$binpath/indexfasta.py $outfile"
    done
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

binpath=$DATADIR3/wk/MPTopo/src/
basename=
mode=
path_predtopo=
path_signalp=


isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        echo Error! Wrong argument: $1 >&2; exit
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -path_predtopo|--path_predtopo) path_predtopo=$2;shift;;
            -path_signalp|--path_signalp) path_signalp=$2;shift;;
            -basename|--basename) basename=$2;shift;;
            -mode|--mode) mode=$2;shift;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        echo Error! Wrong argument: $1 >&2; exit
    fi
    shift
done


if [ "$path_predtopo" == "" ];then
    echo "path_predtopo not set. exit"
    exit 1
fi
if [ "$path_signalp" == "" ];then
    echo "path_signalp not set. exit"
    exit 1
fi
if [ "$basename" == "" ];then
    echo "basename not set. exit"
    exit 1
fi


if [ "$mode" == "" ];then
    case $path_predtopo in 
        *topcons_single*)
            mode=tps
            ;;
        *topcons*)
            mode=tp
            ;;
    esac
fi


signalpfile=$path_signalp/${basename}.signalp_list
if [ ! -f "$signalpfile" ];then
    echo "Warning! $signalpfile does not exist"
fi

case $mode in 
    tps) ProcessPredTopo_tps ;;
    tp) ProcessPredTopo_tp ;;
esac

