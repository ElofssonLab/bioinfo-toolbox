#!/bin/bash

usage="
Usage:  drawReorderedTopoMSA.sh pfamid OR -l pfamidlistfile
    Display multiple alignment of MT protein topologies in png file
Options:
  -outpath  <dir>         set ouput path
  -msapath <dir>          set msa path
  -orderlistpath <dir>    set orderlist file path
  -l       <file>  : set the idListFile
  -text   yes|no   : whether to draw text, default=yes 
  -runtime yes|no  : whether to show running time for each pfamid, default=yes
  -q               : quiet mode
  -h|--help        : print this help message and exit

Created 2012-03-19, updated 2012-03-19, Nanjiang Shu 

Examples:
    drawReorderedTopoMSA.sh PF00122 -msapath msa -orderlistpath msa -outpath msa
"
PrintHelp() {
    echo "$usage"
}

DrawReorderedTopoMSA() { #{{{
    local id=$1
    msafile=$msapath/$id.sorted.orig.topomsa.fa
    reorderedmsafile=$outpath/$id.reordered.topomsa.fa
    orderlistfile=$orderlistpath/${id}-listorder.txt
    if [ ! -s "$msafile" ]; then
        echo "msafile $msafile does not exist or empty. Ignore $id" >&2
        return 1
    fi
    if [ ! -s "$orderlistfile" ]; then
        echo "orderlistfile $orderlistfile does not exist or empty. Ignore $id" >&2
        return 1
    fi
    echo "python $binpath/reordermsa.py -msafile $msafile -orderlist $orderlistfile -o $reorderedmsafile"
    python $binpath/reordermsa.py -msafile $msafile -orderlist $orderlistfile -o $reorderedmsafile

    if [ -f "$reorderedmsafile" ]; then
        python $binpath/drawMSATopo.py -sep n -text n -outpath $outpath\
            $reorderedmsafile
            
        /usr/bin/convert -thumbnail 200 \
            $outpath/$id.reordered.topomsa.png \
            $outpath/thumb.$id.reordered.topomsa.png

        # create full size image with text (text is amino acid sequence)
        python $binpath/drawMSATopo.py -sep n -text y -outpath $outpath\
            $reorderedmsafile
        rm -f $reorderedmsafile
    fi
}

#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

rundir=`dirname $0`
binpath=$rundir

isQuiet=false
outpath=./
msapath=../pfamAna/topoana_mode2_kalignp_nr95
orderlistpath=../pfamAna/figtree_kalignp_nr95
idListFile=
idList=
isDrawText=yes
isShowRunningTime=yes
method_comparison=1

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        idList="$idList $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -msapath|--msapath) msapath=$2;shift;;
            -orderlistpath|--orderlistpath) orderlistpath=$2;shift;;
            -text|--text) isDrawText=$2;shift;;
            -runtime|--runtime) isShowRunningTime=$2;shift;;
            -l|--l|-list|--list) idListFile=$2;shift;;
            -q) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        idList="$idList $1"
    fi
    shift
done


if [ "$idList" == "" -a "$idListFile" == "" ]; then
    echo "$0: Error, input not set!" >&2
    exit
fi


mkdir -p $outpath

res1=
res2=
isShowRunningTime=${isShowRunningTime,,} # get lower case
for id in $idList; do
    if [[ "$isShowRunningTime" == y*  ]]; then 
        res1=$(/bin/date +%s.%N)
    fi
    DrawReorderedTopoMSA $id
    if [[ "$isShowRunningTime" == y*  ]]; then 
        res2=$(/bin/date +%s.%N)
        printf "Running time of %s for drawtopomsa.sh: %8.3F\n" $id  $(echo "$res2 - $res1"|/usr/bin/bc )
    fi
done

if [ "$idListFile" != ""  ]; then
    if [ -f "$idListFile" ] ; then
        for id in $(cat $idListFile); do 
            if [[ "$isShowRunningTime" == y*  ]]; then 
                res1=$(/bin/date +%s.%N)
            fi
            DrawReorderedTopoMSA $id
            if [[ "$isShowRunningTime" == y*  ]]; then 
                res2=$(/bin/date +%s.%N)
                printf "Running time of %s for drawtopomsa.sh: %8.3F\n" $id  $(echo "$res2 - $res1"|/usr/bin/bc )
            fi
        done
    else 
        echo "list file $idListFile does not exist" >&2
    fi
fi

