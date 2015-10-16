#!/bin/bash
# format of the file to make the plot
usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle,-o STR   set the output format (default: eps)
                       can be one of png,term and eps
    -outpath     DIR   set the output path, default: the same as datafile
    -norm              normalize, show relative changes
    -h|--help          print this help message and exit

Created 2012-08-31, updated 2012-09-21, Nanjiang Shu 
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
ylabel="Percentages"
ylabel2="Number of pairs"
isNonOptionArg=false
isNorm=0

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
            -norm|--norm) isNorm=1;;
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

if [ ! -s "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
else
    numDataLine=`awk 'BEGIN{cnt=0}/^[^#]/{cnt++}END{print cnt}' $dataFile`
fi

basename=`basename "$dataFile"` 
scale=(1 1 1 1 1 1 1 1 1)
setyrangeoption1="set yrange[0:100]"
if [ $isNorm -eq 1 ]; then
    ylabel="Relative percentages"
    line1=`awk '{if ($1==0) print}' $dataFile`
    read -ra scale <<< "$line1"
    basename=$basename.norm
    setyrangeoption1="set autoscale y"
fi
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
set ylabel "$ylabel"  font "Arial, 16"
set key top center outside
#set key invert

set style data linespoints
set ytics font 'Times-Roman,16'
set grid noxtics ytics
unset xlabel
set xrange[0:$numDataLine]
$setyrangeoption1
unset xtics

plot \
    '$dataFile' using (\$3/${scale[2]}):xtic(2)   title 'TM2TM' ls 1, \
    ''          using (\$4/${scale[3]}):xtic(2)   title 'TM2GAP' ls 2, \
    ''          using (\$5/${scale[4]}):xtic(2)   title 'TM2SEQ' ls 3

#=============================
#    Plot 2
#=============================
set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.10 
set tmargin at screen 0.30

set xlabel "$xlabel"
set xtics
unset x2tics
set ytics font 'Times-Roman,14'
set ylabel "Number of pairs"
set grid noxtics ytics
set style data histogram
set style fill solid border -1
set xrange[0:$numDataLine]
set  autoscale y
unset key
plot '$dataFile' using 6:xtic(2) ls 2

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

