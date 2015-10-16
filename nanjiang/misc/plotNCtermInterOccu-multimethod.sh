#!/bin/bash
# draw the histogram given the file in the format
# format of the file to make the plot
# #Idx  SeqIDT     Nterm     Cterm  internal   Maximum Occurrence
# 0       0-20    60.919    53.243     3.320    60.919      11778
# 1      20-30    45.620    43.003    13.197    45.620        879
# 2     30-100    58.786    34.238     8.269    58.786        774
usage="
Usage:   plotNCtermInterOccu.sh datafile
Options:
    -outstyle png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2011-11-08, updated 2011-11-08, Nanjiang Shu 
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
            -outstyle|--outstyle) outputStyle=$2;shift;;
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
    echo "Error! dataFile = \"$dataFile\" does not exist or empty. Exit."
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
set terminal png 
set output '$outfile' 
    "
    ;;
    eps)
    outfile=$outpath/$basename.eps
    pdffile=$outpath/$basename.pdf
    outputSetting="
set term postscript eps color
set output '$outfile' 
    "
    ;;
esac

/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"

set title ""
set xlabel "Prediction method"  font "Arial, 16"
set ylabel "Percentages of compared pairs"  font "Arial, 16"
set key top center outside

set style data histogram 
#set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set tics font 'Times-Roman,16'
set xtics scale 0   # remove the bar but keep the tic
#set bmargin 5
set xtics rotate by -30

plot \
    '$dataFile' u 3:xtic(2)   ti 'Unmapped TM at N-terminal' ls 1 , \
    ''          u 4:xtic(2)   ti 'Unmapped TM at C-terminal' ls 2, \
    ''          u 5:xtic(2)   ti 'Unmapped TM at internal regions' ls 3, \
    ''          u (\$1+0.0):(\$6+2):7 ti '' with labels
EOF

echo "plot output to $outfile"


case $outputStyle in
    eps) my_epstopdf $outfile
        echo "Histogram image output to $pdffile"
        if [ -s "$pdffile" ]; then
            rm -f "$outfile"
        fi
        ;;
    *) echo "Histogram image output to $outfile" ;;
esac

