#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
# topcons 78.443  1.502  2.262  0.059  1.896 15.838 68302 
# octopus 79.080  1.630  2.426  0.197  1.813 14.854 66449 
# prodiv 75.885  2.432  4.096  0.505  1.880 15.201 67967 
# pro 78.082  1.630  4.440  0.139  1.851 13.858 61710 
# scampi_msa 79.080  1.630  2.426  0.197  1.813 14.854 66449 
# topcons_single_method4 83.993  3.819  9.461  0.368  2.359  0.000 44511 
# hmmtop 88.234  2.859  6.040  0.528  2.339  0.000 43014 
# memsat 76.714  4.516 14.046  1.863  2.860  0.000 32200 
# scampi_single 85.932  5.505  5.495  0.504  2.564  0.000 39472 
# stmhmm 62.568  4.208 11.262  0.475  2.450 19.036 44158 


usage="
Usage:    plotCmpClass.sh  datafile
Description: 
Options:
    -outstyle|-o png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2011-11-02, updated 2012-08-13, Nanjiang Shu 
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
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

set title ""
set xlabel "Prediction method"  font "Arial, 16"
set ylabel "Percentages of topology comparison classes"  font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
set yrange[0:120]
set label 1 at 3.5, 110
set label 1 "Number of cases"

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
#set bmargin 6
set xtics rotate by -30

plot \
    '$dataFile' using 3:xtic(2)   title 'IDT' ls 1 , \
    ''          using 4:xtic(2)   title 'SHIFT' ls 2, \
    ''          using 5:xtic(2)   title 'INV' ls 3, \
    ''          using 6:xtic(2)   title 'INV_SHIFT' ls 4, \
    ''          using 7:xtic(2)   title 'DUP' ls 5, \
    ''          using 8:xtic(2)   title 'SIGNALP' ls 6, \
    ''          using 9:xtic(2)   title 'Others, same #TM' ls 7,\
    ''          using 10:xtic(2)   title 'Others, different #TM' ls 8,\
    ''          using (\$1):(\$3+\$4+\$5+\$6+\$7+\$8+\$9+\$10+3):11 title '' with labels
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

