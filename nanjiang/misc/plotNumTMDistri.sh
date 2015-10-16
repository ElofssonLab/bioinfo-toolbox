#!/bin/bash
# draw the curves given the file in the format
# #Idx numTM     Nterm     Cterm  internal
# 0        1    17.406    16.867    48.229
# 1        2    15.774    15.465    25.354
# 2        3    14.681    14.597    10.836
# 3        4    12.363    13.470     6.020
# 4        5    11.429    11.519     3.754
# 5        6     7.913     8.151     2.479
# 6        7     6.320     5.564     1.346
# 7        8     4.253     4.206     1.062
# 8        9     3.568     3.512     0.283
# 9       10     1.988     2.905     0.567
# 10      11     2.159     2.255     0.000
# 11      12     1.014     0.679     0.071
# 12      13     0.606     0.434     0.000
# 13      14     0.263     0.087     0.000
# 14      15     0.118     0.130     0.000
# 15      16     0.066     0.014     0.000
# 16      17     0.053     0.072     0.000
# 17      18     0.026     0.043     0.000
# 18      19     0.000     0.029     0.000
# #Total          7595      6919      1412
usage="
Usage:    plotNumTMDistri.sh datafile
Options:
    -outstyle|-o png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2011-11-08, updated 2011-11-09, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=term
outpath=
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

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

if [ ! -s "$dataFile" ]; then 
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
set term postscript eps enhanced color
set output '$outfile' 
    "
    ;;
esac

count[0]=
count[1]=
count[2]=
countstr=`awk '/#Total/ {print $2,$3,$4}' $dataFile`
if [ "$countstr" != "" ];then
    ((i=0))
    for str in $countstr; do 
        count[$i]=$str
        ((i++))
    done
fi

/usr/bin/gnuplot -persist <<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"

set xlabel "Number of TM helices of unmapped regions"  font "Arial, 16"
set ylabel "Normalized frequency of occurrence"  font "Arial, 16"
set key default
set xrange[0:11]

set style data linespoints
set tics font 'Times-Roman,16'
#set bmargin 5

plot \
    '$dataFile' u 3:xtic(2)   ti 'At N-terminal (No. cases: ${count[0]})' ls 1 , \
    ''          u 4:xtic(2)   ti 'At C-terminal (No. cases: ${count[1]})' ls 2, \
    ''          u 5:xtic(2)   ti 'At internal region (No. cases: ${count[2]})' ls 3
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

