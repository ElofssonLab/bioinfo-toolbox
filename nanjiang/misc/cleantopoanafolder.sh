#!/bin/bash 
progname=`basename $0`
usage="
usage:  $progname topoana_folderer

Clean TopoAna folderer

  -h, --help     Print this help message and exit

Created 2013-04-19, updated 2013-04-19, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}
# folderList="
# topoana_clan_mode2_kalignp_topcons
# topoana_clan_mode2_kalignp_topcons_rlty50
# topoana_clan_mode2_kalignp_topcons_rlty75
# topoana_clan_mode2_kalignp_topcons-single
# topoana_pfam_mode2_kalignp_topcons
# topoana_pfam_mode2_kalignp_topcons_rlty50
# topoana_pfam_mode2_kalignp_topcons_rlty75
# topoana_pfam_mode2_kalignp_topcons-single
# "

CleanTopoAnaFolder(){ ##{{{
    local folder="$1"
    local tmpfile=
    tmpfile=$(mktemp /tmp/tmp..clean1.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }
    trap 'rm -f "$tmpfile"' INT TERM EXIT

    # step -1: extract formatted db file 
    local dbIndexfileList=()
    find "$folder" -name "*.indexbin"  > $tmpfile
    if [ -s "$tmpfile" ]; then 
        while read line
        do
            dbIndexfileList+=("$line")
        done < $tmpfile
    fi
    local numFile=${#dbIndexfileList[@]}
    local indexfile=
    local dbname=
    if [ $numFile -gt 0  ]; then
        for ((i=0;i<numFile;i++));do
            indexfile=${dbIndexfileList[$i]}
            dbname=${indexfile%.*}
            $DATADIR3/wk/MPTopo/src/my_extractdb.py  -dbname "$dbname" -splitall -outpath "$folder"
            rm -f "$dbname*.db"
            rm -f "$dbname*.indexbin"
        done
    fi

    #0.
    find "$folder"  -name "*.fa.topo" -print0 | xargs -0 rm -f
    find "$folder"  -name "*.tmp" -print0 | xargs -0 rm -f
    find "$folder"  -name "*.cdhit.nr.fa.tmp*" -print0 | xargs -0 rm -f
    find "$folder"  -name "*.cleaned.topo" -print0 | xargs -0 rm -f
    find "$folder"  -name "*.cleaned.fa" -print0 | xargs -0 rm -f
    find "$folder"  -name "*.homology.cleaned.le1000.fa" -print0 | xargs -0 rm -f
    find "$folder" -type f -name "*.trimmed.topomsa.fa" -print0 | xargs -0 rm -f

    #1.
    find "$folder" -type f -name "*.reordered.topomsa.fa" -print0 | xargs -0 rm -f
    find "$folder" -type l -name "*.reordered.topomsa.fa" -print0 | xargs -0 rm -f
    #2.
    find "$folder" -type f -name "*.cleaned.topomsa.fa" -print0 | xargs -0 rm -f
    find "$folder" -type l -name "*.cleaned.topomsa.fa" -print0 | xargs -0 rm -f
    #3.
    find "$folder" -type f -name "*.renamedid.msa.fasta" -print0 | xargs -0 rm -f
    find "$folder" -type l -name "*.renamedid.msa.fasta" -print0 | xargs -0 rm -f
    #4.
    local file=
    local dir=

    local numfile=
    local rootname=
    find $folder -type f -name "*.clustered.orig.topomsa.fa" > $tmpfile
    if [ -s "$tmpfile" ]; then 
        numfile=`cat $tmpfile | wc -l`
        echo "$numfile are going to be removed."
        for file in $(cat $tmpfile); do 
            rootname=`rootname $file`
            dir=`dirname $file`
            grep "^>" $file > $dir/$rootname.anno
            rm -f $file
        done
    fi
    find "$folder" -type l -name "*.clustered.orig.topomsa.fa" -print0 | xargs -0 rm -f

    #5.
    find "$folder" -type f -name "*.cleaned.topomsa.topowithdgscore" > $tmpfile
    if [ -s "$tmpfile" ]; then 
        numfile=`cat $tmpfile | wc -l`
        echo "$numfile are going to be removed."
        for file in $(cat $tmpfile); do 
            rootname=`rootname $file`
            dir=`dirname $file`
            awk '{if(NR%3== 0 || NR%3 == 1) print}' $file > $dir/$rootname.dgscore
            rm -f $file
        done
    fi
    find "$folder" -type l -name "*.cleaned.topomsa.topowithdgscore" -print0 | xargs -0 rm -f

    case "$folder" in 
        *_clan_*)
            find $folder -name "PF*" -o -name "thumb.PF*" -print0 | xargs -0 rm -f
            ;;
    esac

    #6 gzip multiple sequence alignment
    find "$folder" -type f -name "*.kalignP.fasta" -print0 | xargs -0r gzip -Nf
    find "$folder" -type f -name "*.sorted.orig.topomsa.fa" -print0 | xargs -0r gzip -Nf
    find "$folder" -type f -name "*.topomsa.fa" -print0 | xargs -0r gzip -Nf
    find "$folder" -type f -name "*.origblast.fa" -print0 | xargs -0r gzip -Nf
    rm -f $tmpfile
}
#}}}

folderList=()

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        folderList+=("$1")
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        folderList+=("$1")
    fi
    shift
done

numFolder=${#folderList[@]}
if [ $numFolder -eq 0  ]; then
    echo Input not set! Exit. >&2
    PrintHelp
    exit 1
fi

for ((i=0;i<numFolder;i++));do
    folder=${folderList[$i]}
    CleanTopoAnaFolder "$folder"
done
