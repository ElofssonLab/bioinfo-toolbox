#!/bin/bash

exec_cmd(){
    echo "$*"
    eval "$*"
}

progname=`basename $0`


usage="
Usage: $progname TYPE [TAG]

It should be run within working directory

TYPE:
    all:  Everything
    pro:  Prokaryota
    euk:  Eukaryota
    bac:  Bacteria
    arc:  Archaea

TAG (default: Family):
    Family: For Pfam entries with tag Family
"

TYPE=$1
tag=Family


if [ "$TYPE" == "" ]; then
    echo "TYPE not set"
    echo "$usage"
    exit 1
fi

if [ "$2" != "" ];then
    tag=$2
fi

addname=
case $tag in 
    family|Family)
        addname=.Family
        tag=Family
        ;;
    *)
        addname=""
        ;;
esac

TYPE=`echo $TYPE | tr '[A-Z]' '[a-z]'`
keyword=
addname2=
case $TYPE in 
    al*) keyword=; addname2=;;
    ba*) keyword=Bacteria; addname2=.$keyword;;
    pr*) keyword=Prokaryota; addname2=.$keyword;;
    eu*) keyword=Eukaryota; addname2=.$keyword;;
    ar*) keyword=Archaea; addname2=.$keyword;;
esac


for item in clan pfam; do
    if [ $item == "pfam" -o "$tag"  == "" ];then
        splitpath=split${item}seq_seqfrompfamfasta$addname
        mapfile=../pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20$addname.nr100.filter.fragmented$addname2.${item}id2seqid
        seqdb=../refpro20120604-celluar.selmaxlength-m1.nr100.fasta
        stemname=TMseqdb$addname.nr100.filter.fragmented$addname2
        dbname=$stemname.${item}.fasta
        exec_cmd "$DATADIR3/wk/MPTopo/src/buildFamSeqDB_from_famid2seqidmap.py -mapfile  $mapfile -dbname  $dbname -seqdb $seqdb -tmpdir $splitpath"
        find $splitpath -name "*.fa" | rootname | sort > $stemname.${item}idlist
    fi
done
