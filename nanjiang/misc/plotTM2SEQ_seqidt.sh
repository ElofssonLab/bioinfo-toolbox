#!/bin/bash
# draw histogram of seqidt for TM2SEQ helices
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

Created 2013-06-26, updated 2013-06-26, Nanjiang Shu
"
cmd="$*"
PrintHelp(){
    echo "$usage"
}
MakePlot1(){ #{{{
    local datafile1=$1
    local datafile2=$2
    local b1=$basename
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
    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "purple"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set datafile missing "?"
set style data linespoints 
#set style histogram rowstacked
#set style fill solid border -1
#set boxwidth 0.8
set tics font 'Times-Roman,16'
set autoscale x

plot \
    '$datafile1' using 1:2  title 'Whole sequence' ls 1 , \
    '$datafile2' using 1:2  title 'TM to non-TM segment' ls 2
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
MakePlot2(){ #{{{
    local datafile1=$1
    local datafile2=$2
    local b1=$basename
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
    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "purple"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set datafile missing "?"
#set style fill solid border -1
set style fill transparent solid 0.5 noborder
#set boxwidth 0.8
set tics font 'Times-Roman,16'
set autoscale x

plot \
    '$datafile1' using 1:2  title 'Whole sequence' ls 1 with boxes, \
    '$datafile2' using 1:2  title 'TM to non-TM segment' ls 2 with boxes
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
ylabel="Frequency of occurrences"
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

tmpdatafile1=$(mktemp /tmp/tmp.plotTM2SEQ_seqidt.XXXXXXXXX) || { echo "Failed
to create temp file" >&2; exit 1; }  
tmpdatafile2=$(mktemp /tmp/tmp.plotTM2SEQ_seqidt.XXXXXXXXX) || { echo "Failed
to create temp file" >&2; exit 1; }  
trap 'rm -f "$tmpdatafile1"' INT TERM EXIT
trap 'rm -f "$tmpdatafile2"' INT TERM EXIT

awk '{print $4}' $dataFile | histogram --binsize 3.0 --min 0 > $tmpdatafile1
awk '{print $11}' $dataFile | histogram --binsize 3.0 --min 0 > $tmpdatafile2

cat $tmpdatafile1
cat $tmpdatafile2
#MakePlot1 "$tmpdatafile1" "$tmpdatafile2"
MakePlot2 "$tmpdatafile1" "$tmpdatafile2"
