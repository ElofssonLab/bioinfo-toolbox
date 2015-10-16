#!/bin/bash

# Filename:  add_orgtaxid_to_pdbid.sh
# Description: get the NCBI taxid for the organism by given the pdb id
# Author: Nanjiang Shu (nanjiang.shu@scilifelab.se)

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname [-l LISTFILE] [-o OUTFILE] [PDBID [PDBID...]]
    The input can be four char pdbid or 5 char chain id
    output is 
        pdbid taxid
    tab delimited
OPTIONS:
  -o    OUTFILE  Set output file
  -l    FILE     Set the idListFile, one filename per line
  -h, --help     Print this help message and exit

Created 2014-10-09, updated 2014-10-09, Nanjiang Shu 
"
PrintHelp(){ #{{{
    echo "$usage"
}
#}}}
IsPathExist(){ #{{{
# supply the effective path of the program 
    if ! test -d "$1"; then
        echo Directory \'$1\' does not exist. Aborting $0 >&2
        exit 1
    fi
}
#}}}
GetOrgTaxIDFromPDBID(){ #{{{
    local id=$1
    local pdbfile=`getpdbfilepath $id`
    local txt=`grep "^SOURCE.*ORGANISM_TAXID" $pdbfile | head -n 1`
    local taxid=
    if [ "$txt" != "" ] ; then
        taxid=`echo "$txt" | awk -F "ORGANISM_TAXID: " '{print $2}'  | sed 's/; *$//g'`
        echo -e "$id\t$taxid";
    else
        echo "No organism_taxid info for $id, pdbfile=$pdbfile" >&2
    fi
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=0
outfile=
idListFile=
idList=()

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        idList+=("$1")
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -o|--o) outfile=$2;shift;;
            -l|--l|-list|--list) idListFile=$2;shift;;
            -q|-quiet|--quiet) isQuiet=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        idList+=("$1")
    fi
    shift
done

if [ "$idListFile" != ""  ]; then 
    if [ -s "$idListFile" ]; then 
        while read line
        do
            idList+=("$line")
        done < $idListFile
    else
        echo listfile \'$idListFile\' does not exist or empty. >&2
    fi
fi

numFile=${#idList[@]}
if [ $numFile -eq 0  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi

if [ "$outfile" == "" ];then
    outfile=/dev/stdout
fi

(for ((i=0;i<numFile;i++));do
    id=${idList[$i]}
    GetOrgTaxIDFromPDBID "$id"
done) > $outfile

