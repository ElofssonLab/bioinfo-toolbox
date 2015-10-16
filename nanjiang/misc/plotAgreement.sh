#!/bin/bash
# agreement of different subpredictors
# -1 for missing prediciton
# 0 for non-identical
# 1 for identical
#                 frac       frac     frac
# OCTOPUS     4 -1 0.058    0 0.111   1 0.831
# SCAMPI_seq  5 -1 0.156    0 0.347   1 0.496
# SCAMPI_msa  6 -1 0.058    0 0.111   1 0.831
# PRODIV      7 -1 0.013    0 0.384   1 0.603
# PRO         8 -1 0.324    0 0.243   1 0.433


usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle png|term|eps  set the output, default = terminal
    -outpath DIR            set the output path, default =./
    -h|--help               print this help message and exit

Created 2013-02-21, updated 2013-02-21, Nanjiang Shu
"
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
isNonOptionArg=false
xlabel="Method"
ylabel="Fraction"

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
            -xlabel|--xlabel) xlabel="$2";shift;;
            -ylabel|--ylabel) ylabel="$2";shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFile=$1
    fi
    shift
done

if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
fi

if [ "$outpath" == "" ]; then 
    outpath=`dirname $dataFile`
elif [ ! -d $outpath ]; then
    mkdir -p $outpath
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
    outputSetting="
set term postscript eps color
set output '$outfile' 
    "
    ;;
esac
pdffile="$outpath/$basename.pdf"

/usr/bin/gnuplot -persist<<EOF
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "grey40"
set style line 2 lt 1 lw 1 lc rgb "blue"
set style line 3 lt 1 lw 1 lc rgb "green"

set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set key top center outside
set key invert
set label 1 at 3.5, 110
set label 1 "Number of cases"

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'

plot \
    '$dataFile' using 8:xtic(1)   title 'Identical' ls 3 , \
    ''          using 6:xtic(1)   title 'Non-Identical' ls 2, \
    ''          using 4:xtic(1)   title 'Missing' ls 1
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

