#!/bin/bash

# Filename:   plotCmpClass_mp3_cmpdup_fourfig.sh
# plot A, B, C D, four figure
usage="
Usage: $0 file1 file2 file3 file4
Description: 

    file1: cmpclass for all
    file2: helixcmpclass for all
    file3: cmpclass for Prokaryotes
    file4: cmpclass for eukaryotes

Options:
    -o, -outstyle STR   set the output, (default: eps)
                        can be png, term, eps
    -outpath      DIR   set the output path, default: the same as input file
    -color              Show plots in color, (default: b/w)
    -h, --help          print this help message and exit

Created 2013-10-23, updated 2013-10-23, Nanjiang Shu 

Example:
    $0 file1 file2 file3 file4 -color
"
PrintHelp(){
    echo "$usage"
}

GetYTIC2() { #maxNum
    local maxnum=$1
    #local maxNumExpample=`awk '{print $NF}' $dataFile | sort -rg | head -n 1`
    local ytic2=
    ((ytic2=maxnum/3))
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
MakePlot(){ #{{{  #show also number of cases as multiplot
    local dataFile1="$1"
    local dataFile2="$2"
    local dataFile3="$3"
    local dataFile4="$4"
    local basename=`basename "$dataFile1"` 
    local outputSetting=
    local b1=$basename.fourfig
    local plotOption1=""
    local plotOption2=""

    if [ $isColor -eq 1 ];then
        b1=$b1.color
    fi
    local NF_1=`awk '{print NF}' $dataFile1 | head -n 1` # number of fields
    local NF_2=`awk '{print NF}' $dataFile2 | head -n 1` # number of fields
    local NF_3=`awk '{print NF}' $dataFile3 | head -n 1` # number of fields
    local NF_4=`awk '{print NF}' $dataFile4 | head -n 1` # number of fields
    local maxNumExpample1=`awk '/^[^#]/{print $NF}' $dataFile1 | sort -rg | head -n 1`
    local maxNumExpample2=`awk '/^[^#]/{print $NF}' $dataFile2 | sort -rg | head -n 1`
    local maxNumExpample3=`awk '/^[^#]/{print $NF}' $dataFile3 | sort -rg | head -n 1`
    local maxNumExpample4=`awk '/^[^#]/{print $NF}' $dataFile4 | sort -rg | head -n 1`

    ytic2_1=`GetYTIC2 $maxNumExpample1`
    ytic2_2=`GetYTIC2 $maxNumExpample2`
    ytic2_3=`GetYTIC2 $maxNumExpample3`
    ytic2_4=`GetYTIC2 $maxNumExpample4`

    local N_record_1=`grep "^[^#]" $dataFile1 | wc -l`
    local N_record_2=`grep "^[^#]" $dataFile2 | wc -l`
    local N_record_3=`grep "^[^#]" $dataFile3 | wc -l`
    local N_record_4=`grep "^[^#]" $dataFile4 | wc -l`
    b1=$b1.revert.mtp

    if [ $isColor -eq 0 ];then
        plotOption1="\
        '$dataFile1' using 4:xtic(2)   title 'Inverted'    $fs2 ls 2, \
        ''          using 5:xtic(2)   title 'Duplicated'  $fs3 ls 3, \
        ''          using 6:xtic(2)   title 'TM2GAP'      $fs4 ls 4, \
        ''          using 9:xtic(2)   title 'Mixed'       $fs7 ls 7, \
        ''          using 7:xtic(2)   title 'TM2SEQ'      $fs5 ls 5, \
        ''          using 8:xtic(2)   title 'TM2SP'       $fs6 ls 6  \
        "
        plotOption2="\
        '$dataFile2' using 4:xtic(2)   title 'TM2GAP' $fs2 ls 7, \
        ''          using 5:xtic(2)   title 'TM2SEQ' $fs3 ls 8, \
        ''          using 6:xtic(2)   title 'TM2SP'  $fs4 ls 9  \
        "
        plotOption3="\
        '$dataFile3' using 4:xtic(2)   title 'Inverted'    $fs2 ls 2, \
        ''          using 5:xtic(2)   title 'Duplicated'  $fs3 ls 3, \
        ''          using 6:xtic(2)   title 'TM2GAP'      $fs4 ls 4, \
        ''          using 9:xtic(2)   title 'Mixed'       $fs7 ls 7, \
        ''          using 7:xtic(2)   title 'TM2SEQ'      $fs5 ls 5, \
        ''          using 8:xtic(2)   title 'TM2SP'       $fs6 ls 6  \
        "
        plotOption4="\
        '$dataFile4' using 4:xtic(2)   title 'Inverted'    $fs2 ls 2, \
        ''          using 5:xtic(2)   title 'Duplicated'  $fs3 ls 3, \
        ''          using 6:xtic(2)   title 'TM2GAP'      $fs4 ls 4, \
        ''          using 9:xtic(2)   title 'Mixed'       $fs7 ls 7, \
        ''          using 7:xtic(2)   title 'TM2SEQ'      $fs5 ls 5, \
        ''          using 8:xtic(2)   title 'TM2SP'       $fs6 ls 6  \
        "
    else # in color
        plotOption1="\
        '$dataFile1' using 4:xtic(2)   title 'Inverted'     ls 22, \
        ''           using 5:xtic(2)   title 'Duplicated'   ls 23, \
        ''           using 6:xtic(2)   title 'TM2GAP'       ls 24, \
        ''           using 9:xtic(2)   title 'Mixed'        ls 27, \
        ''           using 7:xtic(2)   title 'TM2SEQ'       ls 25, \
        ''           using 8:xtic(2)   title 'TM2SP'        ls 26  \
        "
        plotOption2="\
        '$dataFile2' using 4:xtic(2)   title 'TM2GAP'   ls 32, \
        ''           using 5:xtic(2)   title 'TM2SEQ'   ls 33, \
        ''           using 6:xtic(2)   title 'TM2SP'    ls 34  \
        "
        plotOption3="\
        '$dataFile3' using 4:xtic(2)   title 'Inverted'     ls 22, \
        ''           using 5:xtic(2)   title 'Duplicated'   ls 23, \
        ''           using 6:xtic(2)   title 'TM2GAP'       ls 24, \
        ''           using 9:xtic(2)   title 'Mixed'        ls 27, \
        ''           using 7:xtic(2)   title 'TM2SEQ'       ls 25, \
        ''           using 8:xtic(2)   title 'TM2SP'        ls 26  \
        "
        plotOption4="\
        '$dataFile4' using 4:xtic(2)   title 'Inverted'     ls 22, \
        ''           using 5:xtic(2)   title 'Duplicated'   ls 23, \
        ''           using 6:xtic(2)   title 'TM2GAP'       ls 24, \
        ''           using 9:xtic(2)   title 'Mixed'        ls 27, \
        ''           using 7:xtic(2)   title 'TM2SEQ'       ls 25, \
        ''           using 8:xtic(2)   title 'TM2SP'        ls 26  \
        "
    fi
    plotOption11="'$dataFile' using $NF_1:xtic(2) ti '' ls 21 "
    plotOption21="'$dataFile' using $NF_2:xtic(2) ti '' ls 21 "
    plotOption31="'$dataFile' using $NF_3:xtic(2) ti '' ls 21 "
    plotOption41="'$dataFile' using $NF_4:xtic(2) ti '' ls 21 "

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
    set term postscript eps color size 7, 6
    set output '$outfile' 
        "
        ;;
    esac

    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
# set style line 1  lt 1 lw 1 lc rgb "red"
# set style line 2  lt 1 lw 1 lc rgb "blue" 
# set style line 3  lt 1 lw 1 lc rgb "green"
# set style line 4  lt 1 lw 1 lc rgb "violet"
# set style line 5  lt 1 lw 1 lc rgb "cyan" 
# set style line 6  lt 1 lw 1 lc rgb "orange"

# color lines
set style line 21  lt 1 lw 1 lc rgb "red"    # IDT
set style line 22  lt 1 lw 1 lc rgb "blue"   # INV
set style line 23  lt 1 lw 1 lc rgb "yellow" # DUP
set style line 24  lt 1 lw 1 lc rgb "green"  # TM2GAP
set style line 25  lt 1 lw 1 lc rgb "violet" # Mixed
set style line 26  lt 1 lw 1 lc rgb "cyan"   # TM2SEQ
set style line 27  lt 1 lw 1 lc rgb "orange" # TM2SP

set style line 32  lt 1 lw 1 lc rgb "navy"   # Helix TM2GAP
set style line 33  lt 1 lw 1 lc rgb "skyblue"    # Helix TM2SEQ
set style line 34  lt 1 lw 1 lc rgb "grey80"  # Helix TM2SP

# black and white
set style line 1  lt 1 lw 1 lc rgb "grey50"
set style line 2  lt 1 lw 1 lc rgb "grey50" 
set style line 3  lt 1 lw 1 lc rgb "grey50"
set style line 4  lt 1 lw 1 lc rgb "grey50"
set style line 5  lt 1 lw 1 lc rgb "grey50" 
set style line 6  lt 1 lw 1 lc rgb "grey50"
set style line 7  lt 1 lw 1 lc rgb "grey50"
set style line 8  lt 1 lw 1 lc rgb "grey50"

set style line 11 lt 1 lw 1 lc rgb "black"
set style line 12 lt 1 lw 1 lc rgb "grey90"


set style line 21 lt 2 lw 2 pt 4 ps 2 lc rgb "black"

set multiplot layout 2, 1 title ""
#=============================
#    Plot 1.1
#=============================

set lmargin at screen 0.08
set rmargin at screen 0.48
set bmargin at screen 0.70
set tmargin at screen 0.95


set title "Protein level" font "Arial, 22"
set label "A)" at screen 0.01, 0.96 font "Arial, 22" 
unset xlabel
unset xtics
set ylabel "$ylabel"  font "Arial, 18"
#set key top center outside
#set key bottom center outside
set key at screen 0.48, 0.94 invert  spacing 1.1
set yrange[0:35]
#set autoscale y

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics font 'Arial,16'
set ytics autofreq
set grid ytics
set autoscale x
set xrange[-1:$N_record_1]

plot $plotOption1

#=============================
#    Plot 1.2
#=============================
set lmargin at screen 0.08
set rmargin at screen 0.48
set bmargin at screen 0.58
set tmargin at screen 0.68

set title ""
set xlabel "$xlabel"  font "Arial, 18"
set xtics font "Arial, 18"
unset x2tics
set ylabel "$ylabel_plot2" font "Arial, 18"
set ytics $ytic2_1 font 'Arial,16'
set grid noxtics ytics
set style data linespoints
#set style data histogram
#set style nofill solid border -1
#set style fill pattern border
set autoscale y
set autoscale y
set yrange[0:]
set autoscale x
set xrange[-1:$N_record_1]
unset key
plot $plotOption11


#=============================
#    Plot 2.1
#=============================

set lmargin at screen 0.59
set rmargin at screen 0.99
set bmargin at screen 0.70
set tmargin at screen 0.95


set title "Helix level" font "Arial, 22"
set label "B)" at screen 0.52, 0.96 font "Arial, 22" 

unset xlabel
#unset ylabel
set ylabel "$ylabel"  font "Arial, 18"
unset xtics
set key at screen 0.985, 0.94 invert spacing 1.1
set yrange[0:7]
#set autoscale y

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics autofreq font 'Arial,16'
set grid ytics
set autoscale x
set xrange[-1:$N_record_2]

plot $plotOption2

#=============================
#    Plot 2.2
#=============================
set lmargin at screen 0.59
set rmargin at screen 0.99
set bmargin at screen 0.58
set tmargin at screen 0.68

set title ""
set xlabel "$xlabel"  font "Arial, 18"
#unset ylabel
set ylabel "$ylabel_plot2" font "Arial, 18"
set xtics font "Arial, 18"
unset x2tics
set ytics font 'Arial,16'
set grid noxtics ytics
set style data linespoints
set autoscale y
set yrange[0:]
set ytics $ytic2_2
set autoscale x
set xrange[-1:$N_record_2]
unset key
plot $plotOption21

#=============================
#    Plot 3.1
#=============================

set lmargin at screen 0.08
set rmargin at screen 0.48
set bmargin at screen 0.20
set tmargin at screen 0.45


set title "Prokaryotes (protein level)" font "Arial, 22"
set label "C)" at screen 0.01, 0.46 font "Arial, 22" 
unset xlabel
unset xtics
set ylabel "$ylabel"  font "Arial, 18"
#set key top center outside
#set key bottom center outside
set key at screen 0.48, 0.44 invert  spacing 1.1
set yrange[0:35]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics autofreq font 'Arial,16'
set grid ytics
set autoscale x
set xrange[-1:$N_record_3]

plot $plotOption3

#=============================
#    Plot 3.2
#=============================
set lmargin at screen 0.08
set rmargin at screen 0.48
set bmargin at screen 0.08
set tmargin at screen 0.18

set title ""
set xlabel "$xlabel"  font "Arial, 18"
set xtics font "Arial, 18"
unset x2tics
set ylabel "$ylabel_plot2" font "Arial, 18"
set ytics font 'Arial,16'
set grid noxtics ytics
set style data linespoints
set autoscale y
set yrange[0:]
set ytics $ytic2_3
set autoscale x
set xrange[-1:$N_record_3]
unset key
plot $plotOption31

#=============================
#    Plot 4.1
#=============================

set lmargin at screen 0.59
set rmargin at screen 0.99
set bmargin at screen 0.20
set tmargin at screen 0.45


set title "Eukaryotes (protein level)" font "Arial, 22"
set label "D)" at screen 0.52, 0.46 font "Arial, 22"
unset xlabel
unset xtics
unset ylabel
set ylabel "$ylabel"  font "Arial, 18"
set key at screen 0.988, 0.44 invert  spacing 1.1
set yrange[0:35]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics autofreq font 'Arial,16'
set grid ytics
set autoscale x
set xrange[-1:$N_record_4]

plot $plotOption4

#=============================
#    Plot 4.2
#=============================
set lmargin at screen 0.59
set rmargin at screen 0.99
set bmargin at screen 0.08
set tmargin at screen 0.18

set title ""
set xlabel "$xlabel"  font "Arial, 18"
unset ylabel
set ylabel "$ylabel_plot2" font "Arial, 18"
set xtics font "Arial, 18"
unset x2tics
set ytics font 'Arial,16'
set grid noxtics ytics
set style data linespoints
set autoscale y
set autoscale y
set yrange[0:]
set ytics $ytic2_4
set autoscale x
set xrange[-1:$N_record_4]
unset key
plot $plotOption41

unset multiplot

EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Histogram image output to $pdffile"
            if [ -s "$pdffile" ]; then
                pdfcrop $pdffile
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
dataFileList=()
xlabel="Sequence identity"
ylabel="Relative frequency (%)"
ylabel_plot2="No. of pairs"
isNonOptionArg=false
isShowNumCase=0
isPlot1=0
isTwoClass=0
isMultiplot=0
title=
isColor=0

fs1="fs solid"
fs2="fs pattern 2"
fs3="fs pattern 5"
fs4="fs pattern 3"
fs5="fs pattern 1"
fs6="fs pattern 4"
fs7="fs pattern 7"

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
            -color|--color) isColor=1;shift;;
            -title|--title) title="$2";shift;;
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

if [ "$numFile" != "4" ]; then 
    exit
fi

dataFile1=${dataFileList[0]}
dataFile2=${dataFileList[1]}
dataFile3=${dataFileList[2]}
dataFile4=${dataFileList[3]}

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile1`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi



stemname=${dataFile1%.*}


MakePlot $dataFile1 $dataFile2 $dataFile3 $dataFile4
