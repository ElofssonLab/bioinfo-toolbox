#!/bin/bash

MakePlot() {
    local dataFile=$1
    local outputSetting=
    local basename=`basename $dataFile`
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
    pdffile=$outpath/$basename.pdf
    ylabel="$ylabel (log scale)"
    /usr/bin/gnuplot -persist<<EOF
$outputSetting
set style line 1 lt 2 lw 1 pt 7 ps 2 lc rgb "black"
set style line 2 lt 3 lw 2 pt 8 ps 2 lc rgb "black"
set style line 3 lt 4 lw 3 pt 9 ps 2 lc rgb "black" 
set style line 4 lt 2 lw 2 pt 7 ps 2 lc rgb "black"
set style line 5 lt 3 lw 1 pt 8 ps 2 lc rgb "black" 

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set xrange[0:13]
set logscale y
#set key top center outside
set datafile missing "?"
set style data linespoints
set tics font 'Times-Roman,16'
plot '$dataFile' using 2:xtic(1) ti 'TM2GAP' ls 1
EOF

    case $outputStyle in
        eps)
            epstopdf $outfile
            echo "Histogram image output to $pdffile"
            ;;
        *)
            echo "Histogram image output to $outfile"
            ;;
    esac
}
#}}}

if [ $# -lt 1 ]; then
    echo "Usage: $0 datafile"
    exit
fi

outputStyle=eps
outpath=
dataFile=
xlabel="Number of TM helices"
ylabel="Occurrences"
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

datadir=`dirname $dataFile`
basename=`basename "$dataFile"`
if [ "$outpath" == "" ]; then
    outpath=$datadir
fi
if [ ! -d "$outpath" ]; then
    mkdir -p "$outpath"
fi

MakePlot $dataFile
