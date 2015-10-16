#!/bin/bash

usage="
Usage: $0 krbiasfile [FILE ...]
Description:
    plot histogram of KRbias distribution, the input is a list of KR-bias values
Options:
    -outstyle|-o png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2013-07-11, updated 2013-07-11, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}
PlotKRBiasHistogram(){
    local dataFile="$1"

    if [ ! -f "$dataFile" ]; then 
        echo "Error! dataFile = \"$dataFile\" does not exist. Ignore."
        return 1
    fi
    local outpath=""

    if [ "$g_outpath" == "" ];then
        outpath=`dirname "$dataFile"`
    else
        outpath=$g_outpath
    fi

    if [ ! -d "$outpath" ];then
        mkdir -p "$outpath"
    fi


    local basename=`basename "$dataFile"` 
    basename=${basename%.*}

    local outputSetting=""
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$outpath/$basename.png
        outputSetting="
set terminal png enhanced
set output '$outfile' 
    "
        ;;
        eps)
        outfile=$outpath/$basename.eps
        pdffile=$outpath/$basename.pdf
        outputSetting="
set term postscript eps enhanced
set output '$outfile' 
    "
        ;;
    esac
    local histfile=$(mktemp /tmp/tmp.plotKRbiasHistogram.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    trap 'rm -f "$histfile"' INT TERM EXIT
    histogram --binsize 1 --min -10.5  -n n -i "$dataFile"  -o "$histfile"
    #histogram --binsize 1  -n n -i "$dataFile"  -o "$histfile"
    maxFreq=`awk 'BEGIN{max=-1} /^[^#]/{if($2>max)max=$2}END{print max}' $histfile`
    ytics=$((maxFreq/5))
    if [ $ytics -le 0 ];then
        ytics=1
    fi

/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lw 1 lc rgb "grey"
set style fill solid  border -1
set title "$title" font "Arial, 28"
set xlabel '$xlabel' font "Arial, 20" 
set ylabel '$ylabel' font "Arial, 20" 
set autoscale y
set tics font "Arial, 20"
set ytics $ytics
set xrange[-10:10]
set xtics 1

set tics font 'Arial,20'
plot '$histfile' using 1:2 title '' ls 1 with boxes

EOF



case $outputStyle in
    eps) my_epstopdf $outfile
        echo "Histogram image output to $pdffile"
        ;;
    *) echo "Histogram image output to $outfile" ;;
esac

}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

xlabel="KR-bias"
ylabel="Occurrences"
outputStyle=eps
g_outpath=
dataFileList=()
title=
isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        dataFileList+=("$1")
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -outstyle|--outstyle|-o) outputStyle=$2;shift;;
            -outpath|--outpath) g_outpath=$2;shift;;
            -title|--title) title="$2";shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFileList+=("$1")
    fi
    shift
done

numFile=${#dataFileList[@]}

if [ $numFile -le 0 ];then
    echo "no datafile set. exit"
    exit 1
fi


for ((i=0;i<numFile;i++));do
    file=${dataFileList[$i]}
    PlotKRBiasHistogram "$file"
done
