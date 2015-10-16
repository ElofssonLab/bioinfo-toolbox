#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
# 0  TM2GAP in Only TM2GAP
# 1  TM2GAP in Both
# 2  TM2SEQ in Only TM2SEQ
# 3  TM2SEQ in Both
# 4  All TM2GAP
# 5  All TM2SEQ
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

Created 2013-04-18, updated 2013-04-18, Nanjiang Shu 
"
cmd="$*"
PrintHelp(){
    echo "$usage"
}
MakePlot() { #{{{ # using lines points, dashed line, black and white, 
    local mode=$1
    local b1=$basename
    local plotOption=
    local tmpdatafile=
    case $mode in 
        only)
            tmpdatafile=$sortedDataFile_only
            plotOption="'$tmpdatafile' using  9 ti 'TM' ls 1"
            plotOption="$plotOption, '' using 10 ti 'nonTM' ls 2"
            b1=$b1.$mode
            ;;
        both)
            tmpdatafile=$sortedDataFile_both
            plotOption="'$tmpdatafile' using 9 ti 'TM' ls 1"
            plotOption="$plotOption, '' using 10 ti 'nonTM' ls 2"
            b1=$b1.$mode
            ;;
        all)
            tmpdatafile=$sortedDataFile_all
            plotOption="'$tmpdatafile' using 9 ti 'TM' ls 1"
            plotOption="$plotOption, '' using 10 ti 'nonTM' ls 2"
            b1=$b1.$mode
            ;;
        *)return 1
            ;;
    esac
    
    #echo "$plotOption"

    local outputSetting=
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
    local pdffile=$outpath/$b1.pdf
    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 pt 2 ps 1 lc rgb "red"
set style line 2 lt 1 lw 1 pt 2 ps 1 lc rgb "blue"
set style line 3 lt 4 lw 3 pt 9 ps 2 lc rgb "black" 
set style line 4 lt 2 lw 2 pt 7 ps 2 lc rgb "black"
set style line 5 lt 3 lw 1 pt 8 ps 2 lc rgb "black" 

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
#set key top center outside
set datafile missing "?"
set style data points
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

}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
xlabel=""
ylabel="DG Value"
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
tmpdir=$(mktemp -d /tmp/tmpdir.plotTM2SEQ_DGvalue.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }  
trap 'rm -rf "$tmpdir"' INT TERM EXIT
sortedDataFile_all=$tmpdir/sorted.all.txt
sortedDataFile_only=$tmpdir/sorted.only.txt
sortedDataFile_both=$tmpdir/sorted.both.txt

sort -k9,9g $dataFile  > $sortedDataFile_all
awk '{if($3=="TM2SEQ")print}' $sortedDataFile_all > $sortedDataFile_only 
awk '{if($3=="TM2GAP_AND_TM2SEQ")print}' $sortedDataFile_all > $sortedDataFile_both

for item in only both all; do
    MakePlot $item
done

rm -rf $tmpdir
echo $tmpdir
