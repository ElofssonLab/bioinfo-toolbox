#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
# #Idx  SeqIDT        IDT     SHIFT       INV INV_SHIFT      DIFF   Maximum Occurrence
# 0       0-10     0.670     2.392     0.315     0.881    95.741    95.741       7608
# 1      10-20     1.881     2.380     0.599     0.932    94.209    94.209       6009
# 2      20-30    41.863    10.515     3.137     1.176    43.309    43.309       4080
# 3      30-40    82.893     5.251     4.987     0.446     6.423    82.893       6056
# 4      40-50    87.557     2.915     5.889     0.211     3.428    87.557       6626
# 5      50-60    90.386     2.434     4.951     0.260     2.008    90.386       7271
# 6      60-70    92.174     2.809     3.836     0.064     1.126    92.174       7795
# 7      70-80    93.070     2.352     3.285     0.062     1.232    93.070       8037
# 8      80-90    92.904     2.942     2.288     0.068     1.798    92.904       7342
# 9     90-100    96.350     1.857     0.944     0.094     0.755    96.350       3178


usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle png|term|eps  set the output, default = terminal
    -outpath DIR            set the output path, default =./
    -h|--help               print this help message and exit

Created 2011-11-02, updated 2011-11-08, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=./
dataFile=
isNonOptionArg=false
xlabel="DG value difference"
ylabel="Frequency"

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

mkdir -p $outpath
#dataFile=/tmp/t2.txt

if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
fi

basename=`basename "$dataFile"` 

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
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

set title "$basename" font "Arial, 26"
set xlabel "$xlabel"  font "Arial, 26"
set ylabel "$ylabel"  font "Arial, 26" 
set nokey
set ytics font 'Times-Roman,26'
set xtics font 'Times-Roman,26' 

set style data histogram 
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,26'

plot '$dataFile' using 2:xtic(1) ti '' ls 1
EOF

case $outputStyle in
    eps) my_epstopdf $outfile
        echo "Histogram image output to $pdffile"
        ;;
    *) echo "Histogram image output to $outfile" ;;
esac

