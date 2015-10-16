#!/bin/bash

usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle png|term|eps  set the output, default = terminal
    -outpath DIR            set the output path, default =./
    -h|--help               print this help message and exit

Created 2011-11-02, updated 2013-11-15, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}


MakePlot(){ #{{{

local dataFile=$1
local outputSetting=
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
set style line 1 lt 1 lw 1 pt 4 ps 1 lc rgb "red"
set style line 2 lt 1 lw 1 pt 7 ps 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "black"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set style data points
set yzeroaxis
set xzeroaxis

set tics font 'Times-Roman,16'

plot \
    '$dataFile' using 1  title 'Cytosolic' ls 1 , \
    ''          using 2  title 'Extracellular' ls 2
#    0 ti '' ls 3 
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
}
#}}}
MakePlot_XY(){ #{{{
    local dataFile1=$1
    local dataFile2=$2
    local xlabel="KR-bias of the protein with N-terminus inside the membrane"
    local ylabel="KR-bias of the protein with N-terminus outside of the membrane"
    local outputSetting=
    local b1=$basename.XY
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$outpath/$b1.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$outpath/$b1.eps
        pdffile=$outpath/$b1.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac
/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 pt 4 ps 1 lc rgb "red"
set style line 2 lt 1 lw 1 pt 7 ps 1 lc rgb "green"
set style line 3 lt 1 lw 3 pt 2 ps 2 lc rgb "grey"
set style line 4 lt 1 lw 3 pt 7 ps 1 lc rgb "black"


set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set style data points
set xrange[-4:16]
set yrange[-12:4]
#set autoscale x
#set autoscale y
set xzeroaxis
set yzeroaxis
set key box
set key height 0.5

set tics font 'Times-Roman,16'

plot \
    '$dataFile1' using 1:2 title 'Single-spanning' ls 3 ,\
    '$dataFile2' using 1:2 title 'Multi-spanning' ls 4
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
}
#}}}
MakePlot_Nin_Nout(){ #{{{
    local dataFile1=$1 #Nin histogram
    local dataFile2=$2 #Nout histogram
    local xlabel="KR-bias"
    local ylabel="Occurrences"
    local outputSetting=
    local b1=$basename.NinNout
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$outpath/$b1.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$outpath/$b1.eps
        pdffile=$outpath/$b1.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac
/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 pt 4 ps 1 lc rgb "red"
set style line 2 lt 1 lw 1 pt 7 ps 1 lc rgb "blue"
set style line 3 lt 1 lw 3 pt 2 ps 2 lc rgb "grey"
set style line 4 lt 1 lw 3 pt 7 ps 1 lc rgb "black"


#set style fill transparent solid 0.5 noborder
set style fill transparent pattern 4 bo
set title ""
set xlabel "$xlabel"  font "Arial, 18"
set ylabel "$ylabel"  font "Arial, 18"
set style data linespoints
set autoscale x
set autoscale y
#set yzeroaxis ls 4
set key top right font "Arial, 18" spacing 1.3

set tics font 'Arial,18'

plot '$dataFile1' using 1:2 title 'N-terminus inside' ls 1 w boxes, \
     '$dataFile2' using 1:2 title 'N-terminus outside' ls 2 w boxes

EOF

case $outputStyle in
    eps) my_epstopdf $outfile
        echo "Histogram image output to $pdffile"
        ;;
    *) echo "Histogram image output to $outfile" ;;
esac
}
#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
isNonOptionArg=false
xlabel=""
ylabel="KR-bias"

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


if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile`
elif [ ! -d "$outpath" ]; then
    mkdir -p "$outpath"
fi


if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
fi

basename=`basename "$dataFile"` 
tmpdatafile=$(mktemp /tmp/tmp.plotINV_krbias.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
tmpdatafile1=$(mktemp /tmp/tmp.plotINV_krbias.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
tmpdatafile2=$(mktemp /tmp/tmp.plotINV_krbias.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
tmpfile_Nin=$(mktemp /tmp/tmp.plotINV_krbias.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
tmpfile_Nout=$(mktemp /tmp/tmp.plotINV_krbias.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
trap 'rm -f "$tmpdatafile"' INT TERM EXIT
trap 'rm -f "$tmpdatafile1"' INT TERM EXIT
trap 'rm -f "$tmpdatafile2"' INT TERM EXIT
trap 'rm -f "$tmpfile_Nin"' INT TERM EXIT
trap 'rm -f "$tmpfile_Nout"' INT TERM EXIT

awk '{if($3=="i"){print $10, $11}else{print $11,$10}}'  $dataFile  > $tmpdatafile
# single spanning
awk '{if($5==1){if($3=="i"){print $10, $11}else{print $11,$10}}}'  $dataFile  > $tmpdatafile1
# multi-spanning
awk '{if($5>1){if($3=="i"){print $10, $11}else{print $11,$10}}}'  $dataFile  > $tmpdatafile2

# Nin
awk '{if($3=="i"){print $10}else{print $11}}'  $dataFile  | histogram --binsize 2 -n no > $tmpfile_Nin
# Nout
awk '{if($3=="i"){print $11}else{print $10}}'  $dataFile  | histogram --binsize 2 -n no > $tmpfile_Nout

MakePlot $tmpdatafile
MakePlot_XY $tmpdatafile1 $tmpdatafile2

echo Nin
cat $tmpfile_Nin

echo Nout
cat $tmpfile_Nout
MakePlot_Nin_Nout $tmpfile_Nin $tmpfile_Nout
