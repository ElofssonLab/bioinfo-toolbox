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

Created 2012-09-11, updated 2012-09-11, Nanjiang Shu 
"
cmd="$*"
PrintHelp(){
    echo "$usage"
}
MakePlot() {
    local mode=$1
    local col1=
    local col2=
    local b1=
    if [ "$mode" == "low" ]; then
        if [ $isEvodist -eq 0 ]; then
            col1=21
            col2=22
        else
            col1=17
            col2=18
        fi
        ti1=$titleLow
        b1=$basename.low
    elif [ "$mode" == "high" ]; then
        if [ $isEvodist -eq 0 ]; then
            col1=23
            col2=24
        else
            col1=19
            col2=20
        fi
        ti1=$titleHigh
        b1=$basename.high
    elif [ "$mode" == "all" ]; then
        if [ $isEvodist -eq 0 ]; then
            col1=25
            col2=26
        else
            col1=21
            col2=22
        fi
        ti1=$titleAll
        b1=$basename.all
    fi
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
set style line 1 lt 1 lw 1 lc rgb "orange"
set style line 2 lt 1 lw 1 lc rgb "grey20" 
set style line 3 lt 1 lw 1 lc rgb "purple"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
#set key top center outside
#set key bottom center outside
#set key invert
#set label 1 at 3.5, 110
#set label 1 "Number of cases"
set datafile missing "?"
set style data linespoints 
#set style histogram rowstacked
#set style fill solid border -1
#set boxwidth 0.8
set tics font 'Times-Roman,16'
set bmargin 5

set xrange[-5:7]

plot \
    '$file1' using $col1:$col2  title 'TM2TM' ls 1 , \
    '$file2' using $col1:$col2  title 'TM2GAP' ls 2, \
    '$file3' using $col1:$col2  title 'TM2SEQ' ls 3
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

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFile=
xlabel="{/Symbol D}G value"
ylabel="Frequency of occurrences"
isNonOptionArg=false
titleLow="SeqIDT < 30%"
titleHigh="SeqIDT >= 30%"
titleAll="All"

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
bsname=${basename%TMMap*hist}

isEvodist=0
case "$dataFile" in 
    *evodist*)
        isEvodist=1
        case "$dataFile" in
            *dgscore_noshift*)
                file1=$dirname/${bsname}TMMap0.evodist.dgscore_noshift.hist
                file2=$dirname/${bsname}TMMap1.evodist.dgscore_noshift.hist
                file3=$dirname/${bsname}TMMap2.evodist.dgscore_noshift.hist
                basename=${bsname}TMMap.evodist.merged.dgscore_noshift.hist
                ;;
            *)
                file1=$dirname/${bsname}TMMap0.evodist.dgscore.hist
                file2=$dirname/${bsname}TMMap1.evodist.dgscore.hist
                file3=$dirname/${bsname}TMMap2.evodist.dgscore.hist
                basename=${bsname}TMMap.evodist.merged.dgscore.hist
                ;;
        esac
        ;;
    *)
        isEvodist=0
        case "$dataFile" in
            *dgscore_noshift*)
                file1=$dirname/${bsname}TMMap0.dgscore_noshift.hist
                file2=$dirname/${bsname}TMMap1.dgscore_noshift.hist
                file3=$dirname/${bsname}TMMap2.dgscore_noshift.hist
                basename=${bsname}TMMap.merged.dgscore_noshift.hist
                ;;
            *)
                file1=$dirname/${bsname}TMMap0.dgscore.hist
                file2=$dirname/${bsname}TMMap1.dgscore.hist
                file3=$dirname/${bsname}TMMap2.dgscore.hist
                basename=${bsname}TMMap.merged.dgscore.hist
                ;;
        esac
esac



MakePlot "low"
MakePlot "high"
MakePlot "all"

