#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
#  #Idx th_TM2TM th_TM2GAP    IDT    INV TM2GAP TM2SEQ  TM2GAP_AND_TM2SEQ
# 0         0.1       0.8 20.919  4.246 11.687 19.351             43.798
# 1         0.2       0.8 20.759  4.218 10.721 19.594             44.708
# 2         0.3       0.8 20.629  4.201  9.549 19.962             45.659
# 3        0.33       0.8 20.629  4.201  9.549 19.962             45.659
# 4         0.4       0.8 20.517  4.177  8.240 20.201             46.865
# 5         0.5       0.8 20.323  4.130  6.971 20.350             48.226
# 6         0.6       0.8 20.072  4.057  5.285 20.503             50.083 

usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle|-o png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2012-08-31, updated 2012-09-11, Nanjiang Shu 
"
PrintHelp(){
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
xlabel="Threshold of fraction of TM2GAP"
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

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages of topology comparison classes"  font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
#set label 1 at 3.5, 110
#set label 1 "Number of cases"

set yrange[0:]
set style data linespoints 
#set style histogram rowstacked
#set style fill solid border -1
#set boxwidth 0.8
set tics font 'Times-Roman,16'
set bmargin 5

plot \
    '$dataFile' using 4:xtic(3)   title 'IDT' ls 1 , \
    ''          using 5:xtic(3)   title 'INV' ls 2, \
    ''          using 6:xtic(3)   title 'Only TM2GAP' ls 3, \
    ''          using 7:xtic(3)   title 'Only TM2SEQ' ls 4, \
    ''          using 8:xtic(3)   title 'Both TM2GAP and TM2SEQ' ls 5
#    ''          using (\$1):(\$3+\$4+\$5+\$6+\$7+3):9 title '' with labels
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

