#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot

usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle|-o png|term|eps  set the output, default = terminal
    -outpath <path>            set the output path, default =./
    -numcase y|n               show/not show number of cases, default: not show
    -h|--help                  print this help message and exit

Created 2013-04-15, updated 2013-04-15, Nanjiang Shu 
"
PrintHelp(){
    echo "$usage"
}
MakePlot_shownumcase(){ #{{{  #showing number of cases
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
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac
    pdffile=$outpath/$basename.pdf
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

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages of topology comparison classes"  font "Arial, 16"
set key top center outside
#set key bottom center outside
set label 1 at 3.5, 110
set label 1 "Number of cases"

set style data histogram 
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
#set bmargin 5

plot \
    '$dataFile' using 3:xtic(2)   title 'IDT' ls 1 , \
    ''          using 4:xtic(2)   title 'INV' ls 2, \
    ''          using (\$1):(\$3+\$4+\$5+\$6+\$7+3):9 title '' with labels
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
MakePlot_no(){ #{{{  #not show number of cases
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
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac
    pdffile=$outpath/$basename.pdf
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

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages"  font "Arial, 16"
set key top center outside
#set key bottom center outside

set style data histogram 
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'

plot \
    '$dataFile' using 3:xtic(2)   title 'Topcons' ls 1 , \
    ''          using 4:xtic(2)   title 'Topcons_single' ls 2
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
xlabel="Sequence identity"
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
                shift;
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

if [ $isShowNumCase -eq 1 ]; then
    MakePlot_shownumcase
else
    MakePlot_no
fi

