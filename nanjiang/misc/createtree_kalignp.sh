#!/bin/bash
usage="
usage:  createtree_kalignp.sh -l pfamidListFile [ID [ID...]]

Description: create tree files for pfamana

    -l       LISTFILE Set list file for PfamID
    -dbname  STR      Database name for sequence
                      (default: /data3/wk/MPTopo/pfamAna/topoana_mode2_kalignp_nr95/homology.cleaned.le1000.kalignP.fasta)
    -datapath DIR     Set datapath, e.g. topoana_mode2_kalignp_nr95
    -outpath  DIR     Set dir for output result, (default: ./result)

Created 2012-03-13, updated 2012-04-13, Nanjiang Shu 

Examples: 
    createtree_kalignp.sh -l t1.pfamidlist -datapath dir -outpath outdir
"


CreateTree() {
    local id=$1
    seqMSAFile=$outpath/$id.kalignp.msa.fa
    topoMSAFile=$outpath/$id.sorted.orig.topomsa.fa
    clusteredTopoMSAAnnoFile=$outpath/$id.clustered.orig.topomsa.anno

    #python $binpath/my_extractdb.py -dbname $dbname $id > $seqMSAFile
    cp -uf $datapath/$id.sorted.orig.topomsa.fa $outpath
    cp  -uf $datapath/$id.homology.cleaned.le1000.kalignP.fasta $seqMSAFile
    grep "^>" $datapath/$id.clustered.orig.topomsa.fa > $clusteredTopoMSAAnnoFile

    if [ ! -f "$seqMSAFile" ]; then 
        echo $id: $seqMSAFile does not exist. Ignore
        return 1
    fi
    treeFile=$outpath/$id.kalignp.fasttree
    if [  ! -s $treeFile ]; then 
        $fasttree_bin/FastTree $seqMSAFile > $treeFile
    fi
    $binpath/sortedTopoMSA2colordef.sh $topoMSAFile > $outpath/$id.cmpclass.colordef.txt
    $binpath/clusteredTopoMSA2colordef.sh $clusteredTopoMSAAnnoFile > $outpath/$id.cluster.colordef.txt
    $binpath/sortedTopoMSA2numTMbardef.sh $topoMSAFile > $outpath/$id.numTMdef.txt
    python $binpath/sortedTopoMSA2inside-outside-colordef.py $topoMSAFile > $outpath/$id.ntermstate.colordef.txt
    python $binpath/itol_pfamtree.py -datapath $outpath -outpath $outpath $id

    #./sortedTopoMSA2taxdef.py $id.sorted.orig.topomsa.fa  > $id.taxdef.txt
    echo $id output
}

progname=`basename $0`
rundir=`dirname $0`

if [ "$1" == "" ] ; then 
    echo "$usage"
    exit
fi
commandline="$0 $*"

binpath=$rundir
fasttree_bin=$rundir/../../../bin

idListFile=
outpath=./
dbname=$rundir/../pfamAna/topoana_mode2_kalignp_nr95/homology.cleaned.le1000.kalignP.fasta 
datapath=$rundir/../pfamAna/topoana_mode2_kalignp_nr95

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        idList="$idList $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) echo "$usage"; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -datapath|--datapath) datapath=$2;shift;;
            -l|--l) idListFile=$2;shift;;
            -dbname|--dbname) dbname=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        idList="$idList $1"
    fi
    shift
done

if [ "$idListFile" != "" ]; then 
    if [ -s "$idListFile" ]; then 
        idList="$idList `/bin/cat $idListFile`"
    else
        echo Warning! idListFile \'$idListFile\' does not exist or empty. >&2
    fi
fi

if [ ${#idList} -eq 0 ]; then
    echo No ID set. Exit. >&2
    exit 1
fi

mkdir -p $outpath
echo Start at    `/bin/date`
general_res1=$(/bin/date +%s.%N)

((cnt=0))
for id in $idList; do 
    res1=$(/bin/date +%s.%N)
    CreateTree $id
    res2=$(/bin/date +%s.%N)
    printf "%s: Running time for %s:    %.3F\n" "$progname" \
        $id  $(echo "$res2 - $res1"|/usr/bin/bc )
    echo
    ((cnt++))
done
echo Finished at `/bin/date`
general_res2=$(/bin/date +%s.%N)
printf "%s: Running time for %d items is: %.3F\n" "$progname" $cnt  \
    $(echo "$general_res2 - $general_res1"|/usr/bin/bc )

