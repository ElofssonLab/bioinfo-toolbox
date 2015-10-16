#!/bin/bash
# draw lines showing the frequency of SHIFT, INV, INV_SHIFT, DUP, SIGNALP,
# DIFF1, DIFF2
# #
# #
# format of the file to make the plot
  #Idx  SeqIDT        OK     SHIFT       INV INV_SHIFT       DUP   SIGNALP     DIFF1     DIFF2   Maximum Occurrence
# 0       0-10     0.066     2.072     0.049     0.658     0.954     0.066    23.166    72.969    72.969      12160
# 1      10-20     5.468     3.784     0.535     0.311     1.725     0.558    21.713    65.906    65.906      26717
# 2      20-30    28.004     6.647     1.867     0.440     3.621     0.793     9.786    48.843    48.843      23197
# 3      30-40    67.401     2.655     2.761     0.122     2.297     1.305     1.675    21.784    67.401      24596
# 4      40-50    82.014     0.926     2.224     0.030     0.448     0.748     0.493    13.116    82.014      26349
# 5      50-60    84.273     0.910     2.126     0.023     0.796     0.932     0.227    10.712    84.273       8794
# 6      60-70    88.218     0.433     1.463     0.000     0.623     0.921     0.569     7.774    88.218       3692
# 7      70-80    91.839     0.496     0.766     0.000     0.180     0.541     0.045     6.132    91.839       2218
# 8      80-90    92.068     0.746     0.881     0.000     0.203     0.678     0.068     5.356    92.068       1475

usage="
Usage:    plotCmpClass_line.sh  datafile
Description: 
Options:
    -outstyle|-o png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2012-08-16, updated 2012-08-16, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=./
dataFile=
xlabel="Sequence identity"
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

mkdir -p $outpath
#dataFile=/tmp/t2.txt

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
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
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
set yrange [0:8]

set style data linespoints
#set style histogram rowstacked
#set style fill solid border -1
#set boxwidth 0.8
set tics font 'Times-Roman,16'
set bmargin 5

plot \
    '$dataFile' using 4:xtic(2)   title 'SHIFT' ls 2, \
    ''          using 5:xtic(2)   title 'INV' ls 3, \
    ''          using 6:xtic(2)   title 'INV_SHIFT' ls 4, \
    ''          using 7:xtic(2)   title 'DUP' ls 5, \
    ''          using 8:xtic(2)   title 'SIGNALP' ls 6
#    ''          using 9:xtic(2)   title 'Others, same #TM' ls 7,\
#    ''          using 10:xtic(2)   title 'Others, different #TM' ls 8,\
#    ''          using (\$1):(\$3+\$4+\$5+\$6+\$7+\$8+\$9+\$10+3):12 title '' with labels
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

