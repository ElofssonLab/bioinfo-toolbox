#!/bin/bash
# format of the file to make the plot
usage="
Usage: $0 datafile1 datafile2
Description: 
    plot helixcmpass 
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
datafileList=()
xlabel="Reliability score threshold"
ylabel="Percentages"
isNonOptionArg=false
isNorm=0

while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        datafileList+=("$1")
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
        datafileList+=("$1")
    fi
    shift
done

numFile=${#datafileList[@]}

if [ $numFile -lt 2 ]; then
    echo "too few input file (<2)" >&2
    exit 1
fi

datafile1=${datafileList[0]}
datafile2=${datafileList[1]}
datadir=`dirname $datafile1`
if [ "$outpath" == "" ]; then
    outpath=$datadir
elif [ ! -d "$outpath"  ]; then
    mkdir -p $outpath
fi

if [ ! -s "$datafile1" ]; then 
    echo "Error! datafile1 \"$datafile1\" does not exist. Exit..."
    exit 1
elif [ ! -s "$datafile2" ]; then 
    echo "Error! datafile2 \"$datafile2\" does not exist. Exit..."
    exit 1
else
    numDataLine1=`awk 'BEGIN{cnt=0}/^[^#]/{cnt++}END{print cnt}' $datafile1`
    numDataLine2=`awk 'BEGIN{cnt=0}/^[^#]/{cnt++}END{print cnt}' $datafile2`
fi

basename=`basename "$datafile1"` 
basename=$basename.two
scale1=(1 1 1 1 1 1 1 1 1)
scale2=(1 1 1 1 1 1 1 1 1)
if [ $isNorm -eq 1 ]; then
    ylabel="Relative percentages"
    line1=`awk '{if ($1==0) print}' $datafile1`
    line2=`awk '{if ($1==0) print}' $datafile2`
    read -ra scale1 <<< "$line1"
    read -ra scale2 <<< "$line2"
    basename=$basename.norm
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
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 

set style line 6 lt 2 lw 1 lc rgb "red"
set style line 7 lt 2 lw 1 lc rgb "green"
set style line 8 lt 2 lw 1 lc rgb "blue" 
set style line 9 lt 2 lw 1 lc rgb "violet"
set style line 10 lt 2 lw 1 lc rgb "cyan" 

set multiplot layout 2, 1 title ""

set lmargin at screen 0.10
set rmargin at screen 0.62
set bmargin at screen 0.55
set tmargin at screen 0.95

set title ""
#set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set key right outside
#set key invert

set style data linespoints
set ytics font 'Times-Roman,16'
set grid noxtics ytics
unset xlabel
set autoscale x
set autoscale y
unset xtics

plot \
    '$datafile1' using (\$3/${scale1[2]}):xtic(2)   title 'IDT, T' ls 1, \
    '$datafile2' using (\$3/${scale2[2]}):xtic(2)   title 'IDT, TS' ls 6


#=============================
#    Plot 2
#=============================
set lmargin at screen 0.10
set rmargin at screen 0.62
set bmargin at screen 0.10 
set tmargin at screen 0.50

set xlabel "$xlabel"
set ylabel "$ylabel"
set xtics
unset x2tics
set key right outside
set ytics font 'Times-Roman,16'
set grid noxtics ytics
set style data linespoints
#set style fill solid border -1
set  autoscale x
set  autoscale y


plot \
    '$datafile1' using (\$4/${scale1[3]}):xtic(2)   title 'INV, T' ls 2, \
    '$datafile1' using (\$5/${scale1[4]}):xtic(2)   title 'Only TM2GAP, T' ls 3, \
    '$datafile1' using (\$6/${scale1[5]}):xtic(2)   title 'Only TM2SEQ, T' ls 4, \
    '$datafile1' using (\$7/${scale1[6]}):xtic(2)   title 'Both TM2GAP and TM2SEQ, T' ls 5, \
    '$datafile2' using (\$4/${scale2[3]}):xtic(2)   title 'INV, TS' ls 7,\
    '$datafile2' using (\$5/${scale2[4]}):xtic(2)   title 'Only TM2GAP, TS' ls 8,\
    '$datafile2' using (\$6/${scale2[5]}):xtic(2)   title 'Only TM2SEQ, TS' ls 9,\
    '$datafile2' using (\$7/${scale2[6]}):xtic(2)   title 'Both TM2GAP and TM2SEQ, TS' ls 10
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

