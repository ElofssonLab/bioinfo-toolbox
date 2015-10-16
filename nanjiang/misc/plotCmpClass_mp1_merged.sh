#!/bin/bash
# draw lines showing the frequency of SHIFT, INV, INV_SHIFT, DUP, SIGNALP,
# DIFF1, DIFF2
# #
# #
# format of the file to make the plot

usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle|-o png|term|eps : set the output, (default: eps)
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2012-08-31, updated 2012-09-11, Nanjiang Shu 
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
xlabel="Reliability score threshold"
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

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
fi

basename=`basename "$dataFile"` 
basename=$basename.line

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
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

set multiplot layout 2, 1 title ""

set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.30
set tmargin at screen 0.85

set title ""
#set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages"  font "Arial, 16"
set key top center outside
#set key bottom center outside
#set key invert

set style data linespoints
set ytics font 'Times-Roman,16'
unset xlabel
unset xtics

plot \
    '$dataFile' using 3:xtic(2)   title 'IDT' ls 1, \
    ''          using 4:xtic(2)   title 'INV' ls 2, \
    ''          using 5:xtic(2)   title 'Only TM2GAP' ls 3, \
    ''          using 6:xtic(2)   title 'Only TM2SEQ' ls 4, \
    ''          using 7:xtic(2)   title 'Both TM2GAP and TM2SEQ' ls 5

#=============================
#    Plot 2
#=============================
set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.10 
set tmargin at screen 0.30

set xlabel "$xlabel"
set xtics
set ylabel "Number of pairs"
set style data histogram
set style fill solid border -1
unset key
plot '$dataFile' using 9:xtic(2) ls 2

unset multiplot
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

