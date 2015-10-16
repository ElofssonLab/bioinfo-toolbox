#!/bin/bash
# format of the file to make the plot
# #Idx    RLTY         0         1         2   Maximum Occurrence
# 0          0     0.000     0.000     0.000     0.000          0
# 1          5     0.000     0.000     0.000     0.000          0
# 2         10     0.000     0.000     0.000     0.000          0
# 3         15     0.000     0.000     0.000     0.000          0
# 4         20     0.000     0.000     0.000     0.000          0
# 5         25     0.000     0.000     0.000     0.000          0
# 6         30     0.000     0.000     0.000     0.000          0
# 7         35     0.000     0.000     0.000     0.000          0
# 8         40     0.000     0.000     0.000     0.000          0
# 9         45     0.000     0.000     0.000     0.000          0
# 10        50     0.000     0.000     0.000     0.000          0
# 11        55    44.444    11.111    44.444    44.444          9
# 12        60    91.346     3.205     5.449    91.346        936
# 13        65    73.613    11.174    15.213    73.613       3714
# 14        70    71.007    12.531    16.462    71.007      17987
# 15        75    73.843    12.104    14.053    73.843      74604
# 16        80    74.807    12.263    12.930    74.807     196388
# 17        85    77.881    10.707    11.412    77.881     350787
# 18        90    83.852     8.681     7.466    83.852     622972
# 19        95    89.890     6.498     3.612    89.890     692061
# 20       100    89.138     6.975     3.887    89.138      57365
# #sum            83.655     8.748     7.598    83.655    2016823
usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle png|term|eps     set the output, (default: eps)
    -outpath <path>            set the output path, default =./
    -norm                      normalize, show relative changes
    -h|--help                  print this help message and exit

Created 2012-11-21, updated 2012-11-21, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}
MakePlot_one_panel() { #{{{
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
    pdffile=$outpath/$basename.pdf

/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 2 lw 2 pt 7  ps 2 lc rgb "black"
set style line 2 lt 3 lw 2 pt 8  ps 2 lc rgb "black" 
set style line 3 lt 4 lw 2 pt 9  ps 2 lc rgb "black"
set style line 4 lt 4 lw 2 pt 10 ps 2 lc rgb "violet"
set style line 5 lt 5 lw 1 pt 11 ps 2 lc rgb "cyan" 
set style line 6 lt 1 lw 1 pt 12 ps 2 lc rgb "orange"
set style line 7 lt 1 lw 1 pt 13 ps 2 lc rgb "grey40"
set style line 8 lt 1 lw 1 pt 14 ps 2 lc rgb "grey70"


set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set key top center #outside

set style data linespoints
set ytics font 'Times-Roman,16'
set grid noxtics ytics
set autoscale x
set autoscale y

plot \
    '$tmpdatafile' using (\$3/${scale[2]}):xtic(2)   title 'TM2TM' ls 1, \
    ''          using (\$4/${scale[3]}):xtic(2)   title 'TM2GAP' ls 2, \
    ''          using (\$5/${scale[4]}):xtic(2)   title 'TM2SEP' ls 3
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

} #}}}
MakePlot_two_panel() { #{{{
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
set key top center #outside
#set key bottom center outside
#set key invert

set style data linespoints
set ytics font 'Times-Roman,16'
set grid noxtics ytics
unset xlabel
set xrange[0:$numDataLine]
$setyrangeoption1
unset xtics

plot \
    '$dataFile' using (\$3/${scale[2]}):xtic(2)   title 'IDT' ls 1, \
    ''          using (\$4/${scale[3]}):xtic(2)   title 'INV' ls 2, \
    ''          using (\$5/${scale[4]}):xtic(2)   title 'Only TM2GAP' ls 3, \
    ''          using (\$6/${scale[5]}):xtic(2)   title 'Only TM2SEQ' ls 4, \
    ''          using (\$7/${scale[6]}):xtic(2)   title 'Both TM2GAP and TM2SEQ' ls 5

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
set ylabel "Number of pairs"
set ytics font 'Times-Roman,14'
set grid noxtics ytics
set style data histogram
set style fill solid border -1
set xrange[0:$numDataLine]
set  autoscale y
unset key
plot '$dataFile' using 8:xtic(2) ls 2

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

} #}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
xlabel="Reliability score"
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
scale=(1 1 1 1 1 1 1 1 1 1)
setyrangeoption1="set yrange[0:100]"
if [ $isNorm -eq 1 ]; then
    ylabel="Relative percentages"
    line1=`awk '{if ($1==0) print}' $dataFile`
    read -ra scale <<< "$line1"
    basename=$basename.norm
    setyrangeoption1="set autoscale y"
fi
basename=$basename.line

tmpdatafile=$(mktemp /tmp/tmp.plotCmpClass_mp1_psbin.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
trap 'rm -f "$tmpdatafile"' INT TERM EXIT

awk '/^[^#]/ {
count=$NF; 
if(count<5000){
    print $1,$2,"? ? ? ? ?"
} else {
    print $0
}
}' $dataFile > $tmpdatafile

#cat $tmpdatafile
MakePlot_one_panel

rm -f $tmpdatafile
