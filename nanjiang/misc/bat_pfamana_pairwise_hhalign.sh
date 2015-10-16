#!/bin/bash

AddAbsolutePath(){ #$fileordir#{{{
    #convert a path to absolute path, 
    #return "" if the file or path does not exist
    if [ ! -e "$1" ]; then 
        return 1
    fi
    curr="$PWD"
    if [ -d "$1" ]; then
        cd "$1" && echo "$PWD"       
    else 
        file="$(basename "$1")"
        cd "$(dirname "$1")" && dir="$PWD"
        echo "$dir/$file"
    fi
    cd "$curr"
    return 0
}
#}}}
exec_cmd(){ #{{{
    echo "$*"
    eval "$*"
} #}}}

#exit
PrintHelp(){
    echo "$usage"
}

usage="
usage: $0 -workdir DIR [-force]
    -force            Force overwrite
    -coverage FLOAT   coverage threshold for hhalign alignment, (default: 0.0)
    -mact     FLOAT   set hhalign mact option, typical is zero
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
    -mp      INT     pairwise comparison method, (default: 3)


    -stopbeforehhalign  Whether stop before running hhalign, just create 
                        pair list

Examples:
    $0 -workdir topcons_agreement_55 -coverage 0 -mact 0
    $0 -workdir topcons_agreement_55 -tax pro -pfamtype Family -coverage 0

Created 2013-06-14, updated 2013-10-09, Nanjiang Shu 
"

rundir=`dirname $0`
rundir=`AddAbsolutePath $rundir`

pfamtype=Family
binpath=$DATADIR3/wk/MPTopo/src
all_seqdb=$DATADIR3/data/uniprot/reference_proteome/refpro20120604.fasta
rootdir=$DATADIR3/wk/MPTopo/pfamAna_refpro/

isNonOptionArg=0
isForceOverWrite=0

MAXPAIR=100
workdir=
maxnum=10000
coverage=0
taxgroup=all
seqtype=all
topotype=all
isStopBeforeHHalign=0
pairwise_comparison_method=3
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -maxnum|--maxnum) maxnum=$2;shift;; 
            -force|--force) isForceOverWrite=1;;
            -maxpair|--maxpair) MAXPAIR=$2;shift;; 
            -workdir|--workdir) workdir=$2;shift;;
            -coverage|--coverage) coverage=$2;shift;;
            -mact|--mact) mact=$2;shift;;
            -tax|--tax) taxgroup=$2;shift;;
            -pfamtype|--pfamtype) pfamtype=$2;shift;;
            -topotype|--topotype) topotype=$2;shift;;
            -seqtype|--seqtype) seqtype=$2;shift;;
            -mp|--mp) pairwise_comparison_method=$2;shift;;
            -stopbeforehhalign|--stopbeforehhalign)isStopBeforeHHalign=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        echo Error! Wrong argument: $1 >&2; exit;
    fi
    shift
done


if [ "$workdir" == "" ]; then
    echo "workdir not set, exit" >&2
    exit 1
fi

if [ "$MAXPAIR" == "" ]; then
    echo "MAXPAIR not set, exit" >&2
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


foldername=`basename $workdir`
predtopodir=$rootdir
bsname_topcons_single=refpro20120604-celluar.selmaxlength-m1.topcons-single
bsname_topcons=refpro20120604-celluar.selmaxlength-m1.topcons.result
bsname_topcons2=query.result
TYPE=
case $foldername in
    topcons_single_agreement*)
        TYPE=TOPCONS_SINGLE
        topo_agreement_str=`echo $foldername | sed 's/topcons_single_agreement_//'`
        topofile=$predtopodir/pred_topcons_single_method4/${bsname_topcons_single}_topcons_single.m1.agree-${topo_agreement_str}${topotype_addname}.topo
        rltyFile=$predtopodir/pred_topcons_single_method4/${bsname_topcons_single}_topcons_single.rlty
        ;;
    topcons_agreement*)
        TYPE=TOPCONS
        topo_agreement_str=`echo $foldername | sed 's/topcons_agreement_//'`
        topofile=$predtopodir/pred_topcons/${bsname_topcons}_TOPCONS.m1.agree-${topo_agreement_str}${topotype_addname}.topo
        rltyFile=$predtopodir/pred_topcons/${bsname_topcons}_TOPCONS.rlty
        ;;
    topcons2_agreement*)
        TYPE=TOPCONS2
        topo_agreement_str=`echo $foldername | sed 's/topcons2_agreement_//'`
        topofile=$predtopodir/pred_topcons/${bsname_topcons2}_TOPCONS_filterSP.m1.agree.m1.agree-${topo_agreement_str}${topotype_addname}.topo
        rltyFile=$predtopodir/pred_topcons/${bsname_topcons2}_TOPCONS.rlty
        ;;
    *)
        TYPE=
        ;;
esac

if [ "$TYPE" == "" ]; then
    echo "wrong type, check workdir name. exit" >&2
    exit 1
elif [ ! -d "$workdir" ] ; then
    mkdir -p "$workdir"
fi

parentdir=$workdir/..
workdir=`AddAbsolutePath $workdir`
parentdir=`AddAbsolutePath $parentdir`

hhalignpath=$workdir/../hhalign
hhprofilepath=$rootdir/hhprofile



signalpfile=$rootdir/pred_signalp/refpro20120604-celluar.selmaxlength-m1.nr100.signalp_list

addname=
addname_pfamtype=
if [ "$pfamtype" != "all" ];then
    addname=${addname}.$pfamtype
    addname_pfamtype=${addname_pfamtype}.$pfamtype
fi

addname=${addname}.nr100.filter.fragmented

if [ "$seqtype" != "all" ];then
    addname=${addname}.$seqtype
fi

if [ "$taxgroup" != "all" ];then
    addname=${addname}.$taxgroup
fi

# prerequisition data files
#=====================================
uniqidlistfile=$rootdir/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20${addname}.uniq.pfam.seqidlist

tm_pfamid2seqidFile=$DATADIR3/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20${addname}.pfamid2seqid

tm_pfamidListFile=$DATADIR3/data/pfam/pfam26.0/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20${addname_pfamtype}.pfamidlist
#======================================



rstpairpath=$workdir/hhalign/rstpair_topcons_cov${coverage}
if [ "$mact" != "" ];then
    rstpairpath=${rstpairpath}_mact$mact
    hhalignpath=${hhalignpath}_mact$mact
fi
if [ ! -d $hhalignpath ]; then
    mkdir -p $hhalignpath
fi
if [ ! -d $rstpairpath ]; then
    mkdir -p  $rstpairpath
fi


echo TYPE=$TYPE
echo foldername=$foldername
echo topo_agreement_str=$topo_agreement_str

echo generating seqdb...

bsname_seqdb=`basename $uniqidlistfile`
rtname_seqdb=${bsname_seqdb%.*}
rtname_seqdb=${rtname_seqdb:0:5} # shorten the name

seqdb=$workdir/${rtname_seqdb}.fasta # changed 2014-05-15, seqdb name is the same as fasta file name due to the recent change in indexfasta.py
seqdbfastafile=${seqdb}
parent_seqdbfastafile=${parentdir}/${rtname_seqdb}.fasta

if [ ! -s ${parent_seqdbfastafile} -o $isForceOverWrite -eq 1 ] ; then
    exec_cmd  "$binpath/my_extractdb.py -l $uniqidlistfile -dbname $all_seqdb  -o ${parent_seqdbfastafile} -q"
    exec_cmd "$binpath/indexfasta.py $parent_seqdbfastafile"
fi

currdir=$PWD
relpath=`python -c "import os, sys; print os.path.relpath(sys.argv[1], sys.argv[2])" $parentdir $workdir  `
cd $workdir
# changed here 2014-05-15, since indexfasta.py has been changed, it is not the
# rootname by the base name is used
exec_cmd "ln -sf $relpath/${rtname_seqdb}.fasta ."
exec_cmd "ln -sf $relpath/${rtname_seqdb}.fasta.index ."
exec_cmd "ln -sf $relpath/${rtname_seqdb}.fasta.indexbin ."
exec_cmd "ln -sf $relpath/${rtname_seqdb}.fasta0.db ."
cd $currdir


echo "$0: Start at    `date`"
res1=$(/bin/date +%s.%N)

pairlistFile=$workdir/${rtname_seqdb}.maxpair${MAXPAIR}.pairlist
pairlistWithPfamIDFile=$workdir/${rtname_seqdb}.maxpair${MAXPAIR}.pairlistwithpfamid


if [ "$TYPE" == "TOPCONS" ]; then
    tmprolistfile=$DATADIR3/wk/MPTopo/pfamAna_refpro/pred_topcons/refpro20120604-celluar.selmaxlength-m1.topcons.result_TOPCONS.m1.agree-${topo_agreement_str}.tmproidlist
elif [ "$TYPE" == "TOPCONS_SINGLE" ]; then
    tmprolistfile=$DATADIR3/wk/MPTopo/pfamAna_refpro/pred_topcons_single_method4/refpro20120604-celluar.selmaxlength-m1.topcons-single_topcons_single.m1.agree-${topo_agreement_str}.tmproidlist
elif [ "$TYPE" == "TOPCONS2" ]; then
    tmprolistfile=$DATADIR3/wk/MPTopo/pfamAna_swissprot/pred_topcons2/${bsname_topcons2}_TOPCONS_filterSP.m1.agree-${topo_agreement_str}.tmproidlist
fi

if [ ! -s "$tmprolistfile" ]; then
    echo tmprolistfile $tmprolistfile does not exist or empty. Exit >&2
    exit
fi

#outname=${rtname_seqdb}.mp${pairwise_comparison_method}.max${maxnum}.maxpair${MAXPAIR}.hhalign
outname=${rtname_seqdb}.mp${pairwise_comparison_method}.cmpdup.max${maxnum}.maxpair${MAXPAIR}.hhalign

echo
echo Step 1: select homologous pairs within each pfam
if [ ! -s "$pairlistWithPfamIDFile" -o $isForceOverWrite -eq 1 ]; then
    exec_cmd  "$binpath/genPairWithinPfam.py -l $tm_pfamidListFile -o $pairlistFile\
        -outwithfamid $pairlistWithPfamIDFile  -mapfile $tm_pfamid2seqidFile \
        -tmprolist $tmprolistfile -m 1 -maxpair $MAXPAIR"
fi

echo

if [ $isStopBeforeHHalign -eq 1 ];then
    exit
fi

echo Step 2: get pairaln and table info from hhalign profiles
pairidlistfile=$rstpairpath/$outname.pairidlistfile
missing_hhr_log=$rstpairpath/$outname.missinghhrlog.txt
missing_pairlist_file=$rstpairpath/$outname.missing.pairlist
hhalignfilelistfile=$rstpairpath/$outname.hhalignfilelist
pairalnfile=$rstpairpath/$outname.pairaln
pairalnstatfile=$rstpairpath/$outname.hhalign.stat
tableinfofile=$rstpairpath/$outname.tableinfo
awk -v dir=$hhalignpath '{print  $1 "_" $2 }' $pairlistWithPfamIDFile > $pairidlistfile

exec_cmd "$binpath/idlist2filelist.py -datapath $hhalignpath -ext .hhr -l $pairidlistfile -o $hhalignfilelistfile 2> $missing_hhr_log"

awk '{print $NF}' $missing_hhr_log | awk -F"_" '{print $1,$2}' > $missing_pairlist_file

#exit 0 # enable this to just create pair list

echo Step 2.1: generate hhalign profiles if not existing
exec_cmd "$binpath/run_hhalign_pairwise.py -l $missing_pairlist_file -hhprofile $hhprofilepath -outpath $hhalignpath"

rm -f $pairidlistfile
rm -f $missing_hhr_log
rm -f $missing_pairlist_file


exec_cmd "$binpath/hhalign2pairaln.py -coverage $coverage -l $hhalignfilelistfile -seqdb $seqdb -o $pairalnfile -os $pairalnstatfile  -ot $tableinfofile"



echo "Step 3: using hhsearch to get duplication file" # added 2014-05-15
dupfile=$rstpairpath/$outname.hhsearch.dup
exec_cmd "$binpath/run_hhsearch_pairwise.py -l $tableinfofile  -hhprofile $hhprofilepath -outpath hhsearch -topofile  $topofile -dupfile $dupfile"

topoalnfile=$rstpairpath/$outname.topoaln.fa
echo Step 4: match topology pairwise alignment
exec_cmd "$binpath/matchMSAtopo -msa $pairalnfile  -topo $topofile  -ignore-badseq no -o $topoalnfile"

echo Step 5: get dg values of alignment
exec_cmd "$binpath/getDGvalueTMOfTopo.sh -topo $topoalnfile -aa $seqdbfastafile -outpath $rstpairpath"

echo Step 6: compare topology
paircmpfile=$rstpairpath/$outname.paircmp
badpaircmpfile=$rstpairpath/$outname.badpaircmp
topoalnwithdgscorefile=$rstpairpath/$outname.topoaln.topowithdgscore
dgscorefile=$rstpairpath/$outname.topoaln.dgscorelist
exec_cmd "$binpath/compareMSATopo.py $topoalnwithdgscorefile -mode 0 -mp $pairwise_comparison_method  -localali $pairalnfile -o $paircmpfile -signalp $signalpfile -cmpsp -obad $badpaircmpfile -cmpdup -dupfile $dupfile -seqidttype 1 -tableinfo $tableinfofile"


echo 
echo Step 7: run pairwise analysis
for alignrange in all full; do
    #exec_cmd "$binpath/ana_paircmp_result.py $paircmpfile -outpath $rstpairpath -tableinfo $tableinfofile -dgscore $dgscorefile -mp 1 -topoaln $topoalnfile -rltyinfo $rltyFile -signalp $signalpfile  -seqaln $pairalnfile  -rmsp -seqidttype 1 -alignrange $alignrange"
    exec_cmd "$binpath/ana_paircmp_result.py $paircmpfile -outpath $rstpairpath -tableinfo $tableinfofile -dgscore $dgscorefile -mp $pairwise_comparison_method -topoaln $topoalnfile -rltyinfo $rltyFile -signalp $signalpfile  -seqaln $pairalnfile  -seqidttype 1 -alignrange $alignrange -dupfile $dupfile"
done

exec_cmd "$binpath/tmp_ana_numTMHeatMap_pfam.py $paircmpfile -outname $outname.tmp_TMHeatMap"

#FINISHED
res2=$(/bin/date +%s.%N)
printf "Running time for %s:  %.3F\n" "$0" $(echo "$res2 - $res1"|/usr/bin/bc )
echo "$0: Finished at `date`"
echo
