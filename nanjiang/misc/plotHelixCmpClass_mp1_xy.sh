#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot

usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle png|term|eps  set the output, default = terminal
    -outpath <path>            set the output path, default =./
    -numcase y|n               show/not show number of cases, default: not show
    -h|--help                  print this help message and exit

Created 2012-12-03, updated 2012-12-03, Nanjiang Shu 
"
PrintHelp(){
    echo "$usage"
}
MakePlot(){ #{{{  #not show number of cases
    local dataFile=$1
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

# set style line 1 lt 2 lw 2 pt 7  ps 1 lc rgb "black"
# set style line 2 lt 3 lw 2 pt 8  ps 1 lc rgb "black" 
# set style line 3 lt 4 lw 2 pt 9  ps 1 lc rgb "black"
# set style line 4 lt 4 lw 2 pt 10 ps 2 lc rgb "violet"
# set style line 5 lt 5 lw 1 pt 11 ps 2 lc rgb "cyan" 
# set style line 6 lt 1 lw 1 pt 12 ps 2 lc rgb "orange"
# set style line 7 lt 1 lw 1 pt 13 ps 2 lc rgb "grey40"
# set style line 8 lt 1 lw 1 pt 14 ps 2 lc rgb "grey70"

set style line 1 lt 1 lw 1 lc rgb "orange"
set style line 2 lt 1 lw 1 lc rgb "grey20"
set style line 3 lt 1 lw 1 lc rgb "purple"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages"  font "Arial, 16"
#set key top center outside
set key left
set grid noxtics ytics

set style data points 
set tics font 'Times-Roman,16'

plot \
    '$dataFile' using 2:3  title 'TM2TM' ls 1 , \
    ''          using 2:4  title 'TM2GAP' ls 2, \
    ''          using 2:5  title 'TM2SEQ' ls 3
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

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
xlabel="Reliability score"
isNonOptionArg=false
isShowNumCase=0

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
            -showcase|--showcase) topt=$2; 
                if [ "${topt:0:1}" == "y" -o  "${topt:0:1}" == "Y" ]; then
                    isShowNumCase=1
                else
                    isShowNumCase=0
                fi
                shift
                ;;
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

MakePlot $dataFile

