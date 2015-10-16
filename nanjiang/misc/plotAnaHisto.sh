#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison vs pid (sequence identity)
# #
# #normalized histogram
# #                       gapless                         local                        global
# pidBins     DIFF SHIFT   INV    OK  sum           DIFF SHIFT   INV    OK  sum           DIFF SHIFT   INV    OK  sum
# <10        0.000 1.000 0.000 0.000    4          0.000 1.000 0.000 0.000    4          1.000 0.000 0.000 0.000    4
# 10-20      0.177 0.822 0.002 0.000  572          0.178 0.822 0.000 0.000  572          0.997 0.002 0.000 0.002  572
# 20-30      0.353 0.218 0.420 0.008  119          0.471 0.218 0.303 0.008  119          0.723 0.017 0.261 0.000  119
# 30-40      0.373 0.000 0.613 0.013   75          0.413 0.013 0.560 0.013   75          0.427 0.000 0.560 0.013   75
# 40-50      0.160 0.187 0.653 0.000   75          0.173 0.187 0.640 0.000   75          0.280 0.120 0.600 0.000   75
# >50        0.052 0.000 0.948 0.000   58          0.052 0.000 0.948 0.000   58          0.052 0.000 0.948 0.000   58
# 

usage="
Usage:     plotAnaHisto.sh datafile
Options:
    -outstyle|-o png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2010-08-18, updated 2010-08-18, Nanjiang
"
function PrintHelp()
{
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=./
dataFile=
isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        dataFile=$1
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -outstyle|--outstyle|-o) outputStyle=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFile=$1
    fi
    shift
done

mkdir -p $outpath

if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
fi

basename=`basename "$dataFile"` 

outputSetting=
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


/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt rgb "black" lw 1    #DIFF
set style line 2 lt rgb "cyan" lw 1     #SHIFT
set style line 3 lt rgb "violet" lw 1     #INV
set style line 4 lt rgb "green" lw 1    #OK

set title "Topology comparison vs sequence identity"
set xlabel "Sequence identity" 0.0,-0.8 # move the xlabel 1 letter downwards
set ylabel "Normalized frequence"
set style data histogram
set key right outside box
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set xtic rotate by -45
#set tics font "Times-Roman,8"
#set xtics font  "Helvetica,8"
#set xtics font 'arial,7'
set tics font 'Verdana,7'
set bmargin 5

plot \
newhistogram "Gapless" , \
'$dataFile' using 2:xtic(1) title col ls 1, '' u 3 ti col ls 2, '' u 4 ti col ls 3, '' u 5 ti col ls 4, \
newhistogram "Local" , \
'$dataFile' using 7:xtic(1) title col ls 1, '' u 8 ti col ls 2, '' u 9 ti col ls 3, '' u 10 ti col ls 4, \
newhistogram "Global" , \
'$dataFile' using 12:xtic(1) title col ls 1, '' u 13 ti col ls 2, '' u 14 ti col ls 3, '' u 15 ti col ls 4
EOF



case $outputStyle in
    eps) my_epstopdf $outfile
        echo "Histogram image output to $pdffile"
        if [ -s "$pdffile" ]; then
            rm -f "$outfile"
        fi
        ;;
    *) echo "Histogram image output to $outfile" ;;
esac

