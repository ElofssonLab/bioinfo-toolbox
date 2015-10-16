#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
#  #Idx  SeqIDT       IDT       INV    TM2GAP    TM2SEQ TM2GAP_AND_TM2SEQ   Maximum Occurrence
# 0       0-10    23.378     4.501    15.654    13.330    43.138    43.138      59702
# 1      10-20    50.891     4.041    19.435    12.237    13.396    50.891      23560
# 2      20-30    70.897     2.489    11.571    10.671     4.372    70.897      25753
# 3      30-40    76.495     2.398     7.739    10.697     2.670    76.495      18762
# 4      40-50    78.959     2.145     4.021    13.445     1.430    78.959      10071
# 5      50-60    83.021     2.340     2.340    11.119     1.181    83.021       4488
# 6      60-70    86.222     1.385     2.450     9.197     0.746    86.222       2816
# 7      70-80    89.414     0.806     1.666     7.523     0.591    89.414       1861
# 8      80-90    89.466     0.964     1.558     7.789     0.223    89.466       1348
# 9     90-100    95.541     0.398     0.318     3.583     0.159    95.541       1256

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
    local outputSetting=
    local outfile=
    local pdffile=
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
    pdffile=$outpath/$basename.pdf
    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 pt 7  ps 1 lc rgb "red"
set style line 2 lt 2 lw 1 pt 8  ps 1 lc rgb "blue" 
set style line 3 lt 3 lw 1 pt 9  ps 1 lc rgb "green"
set style line 4 lt 4 lw 2 pt 10 ps 1 lc rgb "violet"
set style line 5 lt 5 lw 1 pt 11 ps 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 pt 12 ps 1 lc rgb "orange"
set style line 7 lt 1 lw 1 pt 13 ps 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 pt 14 ps 1 lc rgb "grey70"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
#set key top center outside
set key left 
set grid noxtics ytics

set style data points 
set tics font 'Times-Roman,16'

plot \
    '$dataFile' using 2:3  title 'IDT' ls 1 , \
    ''          using 2:4  title 'INV' ls 2, \
    ''          using 2:5  title 'Only TM2GAP' ls 3, \
    ''          using 2:6  title 'Only TM2SEQ' ls 4, \
    ''          using 2:7  title 'Both TM2GAP and TM2SEQ' ls 5
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
ylabel="Percentages"
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

