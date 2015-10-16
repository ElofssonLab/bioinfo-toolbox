#!/bin/bash
# Filename: plotCmpClass_mp3_twofig.sh
# Draw two figure side by side

usage="
Usage: $0 datafile1 datafile2
Description: 
Options:
    -o, -outstyle STR   set the output, (default: eps)
                        can be png, term, eps
    -outpath      DIR   set the output path, default: the same as input file
    -twoclass           show just identical/different
    -h, --help          print this help message and exit

Created 2013-11-13, updated 2013-11-13, Nanjiang Shu 

Example:
    $0 datafile1 datafile2 -twoclass
"
PrintHelp(){
    echo "$usage"
}

GetYTIC2() { #maxNum #{{{
    local maxnum=$1
    #local maxNumExpample=`awk '{print $NF}' $dataFile | sort -rg | head -n 1`
    local ytic2=
    ((ytic2=maxnum/4))
    local numdigit=${#ytic2}
    #echo "maxnum=$maxnum" >&2
    #echo "ytic2=$ytic2" >&2
    #echo "numdigit=$numdigit" >&2
    local a=$ytic2
    for ((i=1;i<numdigit;i++));do
        ((a=a/10))
    done
    for ((i=1;i<numdigit;i++));do
        ((a=a*10))
    done
    ytic2=$a
    if [ $ytic2 -eq 0 ];then
        ytic2=1
    fi
    echo $ytic2
}
#}}}
MakePlot_no_multiplot(){ #{{{  #show also number of cases as multiplot
    local dataFile1="$1"
    local dataFile2="$2"
    local outputSetting=
    local basename=`basename "$dataFile1"`
    local b1=$basename.twofig
    local plotOption1=""
    local plotOption2=""
    local plotOption3=""
    local plotOption4=""
    local NF_1=`awk '{print NF}' $dataFile1 | head -n 1` # number of fields
    local N_record_1=`grep "^[^#]" $dataFile1 | wc -l`
    local NF_2=`awk '{print NF}' $dataFile2 | head -n 1` # number of fields
    local N_record_2=`grep "^[^#]" $dataFile2 | wc -l`
    b1=$b1.mtp


    local maxNumExpample1=`awk '{print $NF}' $dataFile1 | sort -rg | head -n 1`
    local maxNumExpample2=`awk '{print $NF}' $dataFile1 | sort -rg | head -n 1`
    ytic2_1=`GetYTIC2 $maxNumExpample1`
    ytic2_2=`GetYTIC2 $maxNumExpample2`

    if [ $isTwoClass -eq 0 ];then
        plotOption1="\
        '$dataFile1' using 3:xtic(2)   title 'Identical' ls 1, \
        ''          using 4:xtic(2)   title 'Inverted'  ls 2, \
        ''          using 5:xtic(2)   title 'TM2GAP'    ls 3, \
        ''          using 8:xtic(2)   title 'Mixed'     ls 6, \
        ''          using 6:xtic(2)   title 'TM2SEQ'    ls 4, \
        ''          using 7:xtic(2)   title 'TM2SP'     ls 5  \
        "
        plotOption3="\
        '$dataFile2' using 3:xtic(2)   title 'Identical' ls 1, \
        ''          using 4:xtic(2)   title 'Inverted'  ls 2, \
        ''          using 5:xtic(2)   title 'TM2GAP'    ls 3, \
        ''          using 8:xtic(2)   title 'Mixed'     ls 6, \
        ''          using 6:xtic(2)   title 'TM2SEQ'    ls 4, \
        ''          using 7:xtic(2)   title 'TM2SP'     ls 5  \
        "
    else
        plotOption1="\
        '$dataFile1' using (\$4+\$5+\$6+\$7+\$8):xtic(2) title 'Different topology' ls 11 , \
        ''           using 3:xtic(2)   title 'Identical topology' ls 12\
        "
        plotOption3="\
        '$dataFile2' using (\$4+\$5+\$6+\$7+\$8):xtic(2) title 'Different topology' ls 11 , \
        ''           using 3:xtic(2)   title 'Identical topology' ls 12\
        "
        b1=$b1.2cls
    fi
    plotOption2="'$dataFile1' using $NF_1:xtic(2) ti '' ls 21 "
    plotOption4="'$dataFile2' using $NF_2:xtic(2) ti '' ls 21 "

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
    set term postscript eps color size 7, 2.8
    set output '$outfile' 
        "
        ;;
    esac

    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1  lt 1 lw 1 lc rgb "red"
set style line 2  lt 1 lw 1 lc rgb "blue" 
set style line 3  lt 1 lw 1 lc rgb "green"
set style line 4  lt 1 lw 1 lc rgb "violet"
set style line 5  lt 1 lw 1 lc rgb "cyan" 
set style line 6  lt 1 lw 1 lc rgb "orange"
set style line 7  lt 1 lw 1 lc rgb "grey40"
set style line 8  lt 1 lw 1 lc rgb "grey70"
set style line 11 lt 1 lw 1 lc rgb "black"
set style line 12 lt 1 lw 1 lc rgb "grey90"
set style line 21 lt 2 lw 2 pt 4 ps 2 lc rgb "black"

set multiplot layout 2, 2 title ""

#=============================
#    Plot 1.1
#=============================
set lmargin at screen 0.10
set rmargin at screen 0.50
set bmargin at screen 0.35
set tmargin at screen 0.85


set title "$title1"  font "Arial, 20"
unset xlabel
unset xtics
set label "A)" at screen 0.00, 0.96 font "Arial, 22" 
#set xlabel "$xlabel"  font "Arial, 20"
set ylabel "$ylabel"  font "Arial, 17"
#set key top center outside
#set key bottom center outside
#set key top right box
#set key bottom center outside
#set key at screen 0.9, 0.925 font ",14"
set key at screen 0.23, 0.98 font ",16"
set key spacing 1.15
set key invert

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics font 'Arial,16'
set grid ytics
set autoscale x
set autoscale y
set xrange[-1:$N_record_1]
set yrange[0:100]

plot $plotOption1

#=============================
#    Plot 1.2
#=============================
set lmargin at screen 0.10
set rmargin at screen 0.50
set bmargin at screen 0.15
set tmargin at screen 0.30

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set xtics font "Arial,16"
unset x2tics
set ylabel "$ylabel_plot2" font "Arial, 17"
set ytics font 'Arial,14'
set grid noxtics ytics
set style data linespoints
set autoscale x
set autoscale y
set yrange[0:]
set ytics $ytic2_1
set xrange[-1:$N_record_1]
unset key

plot $plotOption2

#=============================
#    Plot 2.1
#=============================
set lmargin at screen 0.58
set rmargin at screen 0.98
set bmargin at screen 0.35
set tmargin at screen 0.85


set title "$title2"  font "Arial, 20"
set label "B)" at screen 0.53, 0.96 font "Arial, 22" 
unset xlabel
unset ylabel
unset xtics
unset key

set style data histogram
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics font 'Arial,16'
set grid ytics
set autoscale x
set autoscale y
set xrange[-1:$N_record_2]
set yrange[0:100]
set ytic 20

plot $plotOption3

#=============================
#    Plot 2.2
#=============================
set lmargin at screen 0.58
set rmargin at screen 0.98
set bmargin at screen 0.15
set tmargin at screen 0.30

set title ""
set xlabel "$xlabel"  font "Arial, 17"
set xtics font "Arial,16"
unset x2tics
#set ylabel "$ylabel_plot2"
unset ylabel
set ytics font 'Arial,14'
set grid noxtics ytics
set style data linespoints
set autoscale x
set autoscale y
set yrange[0:]
set ytics $ytic2_2
set xrange[-1:$N_record_2]
unset key

plot $plotOption4



unset multiplot

EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            pdfcrop $pdffile
            echo "Histogram image output to $pdffile"
            ;;
        *) echo "Histogram image output to $outfile" ;;
    esac

}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

title1=
title2=
outputStyle=eps
outpath=
dataFileList=()
xlabel="Number of agreed topology predictors"
ylabel="Relative frequency (%)"
ylabel_plot2="No. of pairs"
isNonOptionArg=false
isShowNumCase=0
isPlot1=0
isTwoClass=0
isMultiplot=0

fs1="fs solid"
fs2="fs pattern 2"
fs3="fs pattern 3"
fs4="fs pattern 1"
fs5="fs pattern 4"
fs6="fs pattern 7"

while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        dataFileList+=($1)
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -outstyle|--outstyle|-o) outputStyle=$2;shift;;
            -xlabel|--xlabel) xlabel="$2";shift;;
            -ylabel|--ylabel) ylabel="$2";shift;;
            -title1|--title1) title1="$2";shift;;
            -title2|--title2) title2="$2";shift;;
            -plot1|--plot1) isPlot1=1;;
            -twoclass|--twoclass) isTwoClass=1;;
            -multiplot|--multiplot) isMultiplot=1;;
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
        dataFileList+=($1)
    fi
    shift
done

numFile=${#dataFileList[@]}
if [ "$numFile" != "2" ]; then 
    echo "Error! numFile ($numFile) != 2. Exit..."
    exit
fi
dataFile1=${dataFileList[0]}
dataFile2=${dataFileList[1]}

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile1`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi



stemname=${dataFile1%.*}
stemname=${stemname}.twofig


MakePlot_no_multiplot "$dataFile1" "$dataFile2"

