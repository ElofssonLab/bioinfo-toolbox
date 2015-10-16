#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle,-o STR    Set the output, (default: eps)
                        can be png, term, eps
    -outpath DIR        Set the output path, (default: ./)
    -xlabel STR         Set xlabel, (default: Sequence identity)
    -ylabel STR         Set Ylabel
    -h,--help           Print this help message and exit

Created 2013-06-25, updated 2013-06-25, Nanjiang Shu 
"
cmd="$*"
PrintHelp(){
    echo "$usage"
}
MakePlot1(){ #{{{ # using histogram
    local dataFile=$1
    local plotOption=""
    local b1=$basename

    plotOption="\
        '$dataFile' using 2:xtic(1)  title 'Same side' fs solid 1 ls 1,\
        ''          using 3:xtic(1)  title 'Different side' fs solid 1 ls 2
    "
    outputSetting=
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
            set term postscript eps enhanced
            set output '$outfile' 
            "
            ;;
    esac
    pdffile=$outpath/$b1.pdf
#echo "plot $plotOption"
    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue"
set style line 3 lt 1 lw 1 lc rgb "black" 
set style line 4 lt 1 lw 1 lc rgb "grey70"
set style line 5 lt 1 lw 1 lc rgb "black" 
set style line 6 lt 1 lw 1 lc rgb "grey70"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
#set key top center outside
set datafile missing "?"
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set tics font 'Times-Roman,16'

plot $plotOption
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


rm -f $tmpdatafile
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
xlabel="Difference in the number of N-terminal TM helices"
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


if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi
if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    echo "$cmd"
    exit
fi


basename=`basename "$dataFile"` 
dirname=`dirname $dataFile`

tmpdatafile=$(mktemp /tmp/tmp.plotNumTMHeatMap.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
trap 'rm -f "$tmpdatafile"' INT TERM EXIT

case "$dataFile" in 
    *DiffNumTMNterm*) xlabel="Difference in the number of N-terminal TM helices";;
    *DiffNumTMCterm*) xlabel="Difference in the number of C-terminal TM helices";;
esac

#awk '/^[^#]/{if ($1!=0)print $1-1, $1, $2, $3}' $dataFile > $tmpdatafile
awk '{if ($1!=0)print $1, $2, -$3}' $dataFile > $tmpdatafile

MakePlot1 $tmpdatafile

