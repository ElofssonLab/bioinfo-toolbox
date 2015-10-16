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
    -outpath DIR            set the output path, default: dirname(infile)
    -h|--help               print this help message and exit

Created 2011-11-02, updated 2011-11-08, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}
MakePlot(){ #{{{
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
    pdffile=$outpath/$basename.pdf
/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "grey10"
set style line 2 lt 1 lw 1 lc rgb "grey30"
set style line 3 lt 1 lw 1 lc rgb "grey50" 
set style line 4 lt 1 lw 1 lc rgb "grey70"
set style line 5 lt 1 lw 1 lc rgb "grey90" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set key top center outside
set datafile missing '?'
set autoscale x
set autoscale y
set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.9
set tics font 'Times-Roman,16'

plot \
    '$dataFile' using 3:xtic(2)   title 'SeqIDT 0 - 10%' ls 1 , \
    ''          using 4:xtic(2)   title 'SeqIDT 10 - 20%' ls 2, \
    ''          using 5:xtic(2)   title 'SeqIDT 20 - 30%' ls 3, \
    ''          using 6:xtic(2)   title 'SeqIDT 30 - 80%' ls 4, \
    ''          using 7:xtic(2)   title 'SeqIDT 80 - 100%' ls 5
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
isNonOptionArg=false
xlabel="Reliability score"
ylabel="Number of pairs"

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
    mkdir -p $outpath
fi

if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
fi

basename=`basename "$dataFile"` 

tmpdatafile=$(mktemp /tmp/tmp.plotCmpClass_mp1_psbin.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
trap 'rm -f "$tmpdatafile"' INT TERM EXIT

awk '/^[^#]/ {
count=$NF; 
if(count<1000){
    #print $1,$2,"? ? ? ? ? ?"
} else {
    print $0
}
}' $dataFile > $tmpdatafile

#cat $tmpdatafile
MakePlot $tmpdatafile

rm -f $tmpdatafile

