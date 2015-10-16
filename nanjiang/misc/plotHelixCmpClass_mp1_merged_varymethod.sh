#!/bin/bash
# format of the file to make the plot
# #Idx          method    IDT    INV TM2GAP TM2SEQ TM2GAP_AND_TM2SEQ      Count
# 0            topcons 35.139  4.393 15.285 13.843 31.340      88714
# 1         scampi_msa 35.072  4.749 15.113 13.578 31.488      85904
# 2            octopus 35.072  4.749 15.113 13.578 31.488      85904
# 3        polyphobius 35.162  5.499 15.308 15.300 28.732      86002
# 4             prodiv 36.498  3.917 17.322 12.037 30.225      86023
# 5                pro 34.359  4.068 14.079 16.789 30.704      85069
# 6     topcons_single 25.529  6.675 10.888 22.673 34.235      88357
# 7      scampi_single 23.319  6.686  8.415 23.481 38.099      88524
# 8            phobius 25.436  7.175 10.349 24.541 32.499      82020
# 9             hmmtop 23.307  7.244 11.566 22.099 35.785      87701
# 10            memsat 16.033  8.004  8.925 29.040 37.998      85806
# 11            stmhmm 23.516  5.611  8.925 27.140 34.808      87090
usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle|-o png|term|eps : set the output, (default: eps)
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2012-09-21, updated 2012-09-21, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

PlotLinePoints(){ #{{{
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
#set key invert

set rmargin at screen 0.9
set style data linespoints
set ytics font 'Times-Roman,16'
set grid noxtics ytics
set yrange[0:100]
set xtic rotate by -30
plot \
    '$dataFile' using 3:xtic(2)   title 'TM2TM' ls 1, \
    ''          using 4:xtic(2)   title 'TM2GAP' ls 2, \
    ''          using 5:xtic(2)   title 'TM2SEQ' ls 3
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
PlotHistogram(){ #{{{
basename=`basename "$dataFile"` 
basename=$basename


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


set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages"  font "Arial, 16"
set key top center outside
set key invert
set xtic rotate by -30

set style data histogram
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics font 'Times-Roman,16'
set grid noxtics ytics

plot \
    '$dataFile' using 3:xtic(2)   title 'TM2TM' ls 1, \
    ''          using 4:xtic(2)   title 'TM2GAP' ls 2, \
    ''          using 5:xtic(2)   title 'TM2SEQ' ls 3
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

outputStyle=eps
outpath=
dataFile=
xlabel="Method"
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

PlotLinePoints
PlotHistogram
