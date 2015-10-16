#!/bin/bash
# format of the file to make the plot
usage="
Usage: $0 datafile1 datafile2
Description: 
Options:
    -outstyle|-o png|term|eps   set the output, (default: eps)
    -outpath <path>             set the output path, default =./
    -norm                       normalize, show relative changes
    -h|--help                   print this help message and exit

Created 2012-08-31, updated 2012-09-20, Nanjiang Shu 
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
xlabel="Fraction of number of pairs"
ylabel="Percentages"
ylabel2="Number of pairs"
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

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

numFile=${#datafileList[@]}
if [ $numFile -lt 2 ]; then
    echo "two few input datafiles (<2)" >&2
    exit 1
fi

datafile1=${datafileList[0]}
datafile2=${datafileList[1]}
if [ ! -s "$datafile1" -o ! -s "$datafile2" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
else
    numDataLine1=`awk 'BEGIN{cnt=0}/^[^#]/{cnt++}END{print cnt}' $datafile1`
    numDataLine2=`awk 'BEGIN{cnt=0}/^[^#]/{cnt++}END{print cnt}' $datafile2`
fi

basename=`basename "$datafile1"` 
basename=$basename.two
scale1=(1 1 1 1 1 1 1 1 1)
scale2=(1 1 1 1 1 1 1 1 1)
setyrangeoption1="set autoscale y"
if [ $isNorm -eq 1 ]; then
    ylabel="Relative percentages"
    line1=`awk '{if ($1==0) print}' $datafile1`
    line2=`awk '{if ($1==0) print}' $datafile2`
    read -ra scale1 <<< "$line1"
    read -ra scale2 <<< "$line2"
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
set style line  1 lt 1 lw 1 lc rgb "red"
set style line  2 lt 1 lw 1 lc rgb "blue" 
set style line  3 lt 1 lw 1 lc rgb "green"
set style line  4 lt 1 lw 1 lc rgb "violet"
set style line  5 lt 1 lw 1 lc rgb "cyan" 

set style line  6 lt 2 lw 1 lc rgb "red"
set style line  7 lt 2 lw 1 lc rgb "blue" 
set style line  8 lt 2 lw 1 lc rgb "green"
set style line  9 lt 2 lw 1 lc rgb "violet"
set style line 10 lt 2 lw 1 lc rgb "cyan" 


set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set key top center outside

set style data linespoints
set ytics font 'Times-Roman,16'
set grid noxtics ytics
set xrange[*:*] reverse
set autoscale y

plot \
    '$datafile1' using (\$8/${scale1[7]}):(\$3/${scale1[2]}) title 'IDT, Topcons' ls 1, \
    '$datafile1' using (\$8/${scale1[7]}):(\$4/${scale1[3]}) title 'INV, Topcons' ls 2, \
    '$datafile1' using (\$8/${scale1[7]}):(\$5/${scale1[4]}) title 'Only TM2GAP, Topcons' ls 3, \
    '$datafile1' using (\$8/${scale1[7]}):(\$6/${scale1[5]}) title 'Only TM2SEQ, Topcons' ls 4, \
    '$datafile1' using (\$8/${scale1[7]}):(\$7/${scale1[6]}) title 'Both TM2GAP and TM2SEQ, Topcons' ls 5, \
    '$datafile2' using (\$8/${scale2[7]}):(\$3/${scale2[2]}) title 'IDT, Topcons_single' ls 6, \
    '$datafile2' using (\$8/${scale2[7]}):(\$4/${scale2[3]}) title 'INV, Topcons_single' ls 7, \
    '$datafile2' using (\$8/${scale2[7]}):(\$5/${scale2[4]}) title 'Only TM2GAP, Topcons_single' ls 8,\
    '$datafile2' using (\$8/${scale2[7]}):(\$6/${scale2[5]}) title 'Only TM2SEQ, Topcons_single' ls 9,\
    '$datafile2' using (\$8/${scale2[7]}):(\$7/${scale2[6]}) title 'Both TM2GAP and TM2SEQ, Topcons_single' ls 10
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

