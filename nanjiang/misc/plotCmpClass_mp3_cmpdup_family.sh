#!/bin/bash

# Filename:   plotCmpClass_mp3_cmpdup_family.sh
usage="
Usage: $0 datafile
DESCRIPTION: 
    Plot cmpclass for all families, bin with family
    sorted in descending order by the number of pairs of the family
OPTIONS:
    -o, -outstyle STR   set the output, (default: eps)
                        can be png, term, eps
    -outpath      DIR   set the output path, default: the same as input file
    -color              Show plots in color, (default: b/w)
    -h, --help          print this help message and exit

Created 2015-03-12, updated 2015-03-12, Nanjiang Shu

Example:
    $0 datafile 
"
PrintHelp(){
    echo "$usage"
}
exec_cmd(){
    echo "$*"
    eval "$*"
}

GetYTIC2() { #maxNum#{{{
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
#}}}
CreateDataFile1(){ #{{{
    local infile="$1"
    local outfile="$2"
    exec_cmd "cat $infile | grep -v \"^#\" | sort -k3,3rg $infile | nl > $outfile"
}
#}}}


Makeplot_mtp(){ #{{{  #show also number of cases as multiplot
    local dataFile="$1"
    local isNormalize="$2"
    local isTwoClass="$3"

    local basename=`basename "$dataFile"` 
    local outputSetting=
    local b1=$basename
    local plotOption1=""
    local plotOption2=""
    local NF=`awk '{print NF}' $dataFile | head -n 1` # number of fields
    local N_record=`grep "^[^#]" $dataFile | wc -l`
    b1=$b1.mtp
    if [ "$isNormalize" == "1" ];then
        b1=$b1.norm
    fi


    local maxNumExpample=`awk '/^[^#]/{print $4}' $dataFile | sort -rg | head -n 1`
    ytic2=`GetYTIC2 $maxNumExpample`
    maxNumExpample=`awk '/^[^#]/{print $6}' $dataFile | sort -rg | head -n 1`
    ytic3=`GetYTIC2 $maxNumExpample`
    maxNumExpample=`awk '/^[^#]/{print $5}' $dataFile | sort -rg | head -n 1`
    ytic4=`GetYTIC2 $maxNumExpample`

    if [ $isTwoClass ==  "0" ];then
        if [ "$isNormalize" == "0" ];then
            plotOption1="\
            '$dataFile' using  9:xtic(1)   title 'Identical'    ls 21, \
            ''          using 10:xtic(1)   title 'Inverted'     ls 22, \
            ''          using 11:xtic(1)   title 'Duplicated'   ls 23, \
            ''          using 12:xtic(1)   title 'TM2GAP'       ls 24, \
            ''          using 15:xtic(1)   title 'Mixed'        ls 27, \
            ''          using 13:xtic(1)   title 'TM2SEQ'       ls 25, \
            ''          using 14:xtic(1)   title 'TM2SP'        ls 26\
            "
        else
            plotOption1="\
            '$dataFile' using  (\$9/\$4):xtic(1)   title 'Identical'    ls 21, \
            ''          using (\$10/\$4):xtic(1)   title 'Inverted'     ls 22, \
            ''          using (\$11/\$4):xtic(1)   title 'Duplicated'   ls 23, \
            ''          using (\$12/\$4):xtic(1)   title 'TM2GAP'       ls 24, \
            ''          using (\$15/\$4):xtic(1)   title 'Mixed'        ls 27, \
            ''          using (\$13/\$4):xtic(1)   title 'TM2SEQ'       ls 25, \
            ''          using (\$14/\$4):xtic(1)   title 'TM2SP'        ls 26\
            "
        fi
    else
        if [ "$isNormalize" == "0" ];then
            plotOption1="\
            '$dataFile'  using (\$9/\$4):xtic(2)   title 'Identical topology' ls 12,\
            ''           using (\$10+\$11+\$12+\$13+\$14+\$15):xtic(1) title 'Different topology' ls 11\
            "
        else
            plotOption1="\
            '$dataFile'   using (\$9/\$4):xtic(2)   title 'Identical topology' ls 12,\
            ''            using ((\$10+\$11+\$12+\$13+\$14+\$15)/\$4):xtic(1) title 'Different topology' ls 11\
            "
        fi
        b1=$b1.2cls
    fi
    #plotOption2="'$dataFile' using $3:xtic(1) ti '' ls 31 "
    plotOption2="'$dataFile' using 1:4 ti '' ls 31 " # numPair
    plotOption3="'$dataFile' using 1:6 ti '' ls 31 " # numSeq
    plotOption4="'$dataFile' using 1:5 ti '' ls 31 " # numPredicted TMprotein

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

    /usr/bin/gnuplot -persist<<EOF 
$outputSetting

# color line style
set style line 21  lt 1 lw 1 lc rgb "red"    # IDT
set style line 22  lt 1 lw 1 lc rgb "blue"   # INV
set style line 23  lt 1 lw 1 lc rgb "yellow" # DUP
set style line 24  lt 1 lw 1 lc rgb "green"  # TM2GAP
set style line 25  lt 1 lw 1 lc rgb "violet" # Mixed
set style line 26  lt 1 lw 1 lc rgb "cyan"   # TM2SEQ
set style line 27  lt 1 lw 1 lc rgb "orange" # TM2SP


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


set style line 31 lt 1 lw 2 pt 6 ps 0 lc rgb "blue"



set multiplot layout 4, 1 title ""

set lmargin at screen 0.15
set rmargin at screen 0.80
set bmargin at screen 0.50
set tmargin at screen 0.95


set title "$title"
unset xlabel
unset xtics
#set xlabel "$xlabel"  font "Arial, 14"
set ylabel "$ylabel"  font "Arial, 14"
#set key top center outside
#set key bottom center outside
#set key at screen 0.9, 0.925
set key at screen 0.975, 0.95
#set key top left outside
set key invert
set autoscale y
#set yrange[0:100]

set style data histogram 
set style histogram rowstacked
set style fill solid noborder
#set boxwidth 0.8
set ytics font 'Times-Roman, 14'
#set grid ytics
set autoscale x
set xrange[-1:$N_record]

plot \
    $plotOption1

#=============================
#    Plot 2
#=============================
set lmargin at screen 0.15
set rmargin at screen 0.80
set bmargin at screen 0.37
set tmargin at screen 0.47

set title ""
set xlabel "$xlabel"  font "Arial, 14"
unset xlabel
unset xtics
unset x2tics
set ylabel "$ylabel_plot2" font "Arial, 12"  rotate by 0 offset 0
set ytics font 'Times-Roman, 14'
set grid noxtics ytics
set style data linespoints
set autoscale x
set autoscale y
set yrange[0:]
set ytics $ytic2
set xrange[-1:$N_record]
unset key

plot \
    $plotOption2

#=============================
#    Plot 3
#=============================
set lmargin at screen 0.15
set rmargin at screen 0.80
set bmargin at screen 0.24
set tmargin at screen 0.34

set title ""
set xlabel "$xlabel"  font "Arial, 14"
unset xlabel
unset xtics
unset x2tics
set ylabel "$ylabel_plot3" font "Arial, 12"  rotate by 0 offset 1.5,0.75
set ytics font 'Times-Roman, 14'
set grid noxtics ytics
set style data linespoints
set autoscale x
set autoscale y
set yrange[0:]
set ytics $ytic3
set xrange[-1:$N_record]
unset key

plot \
    $plotOption3

#=============================
#    Plot 4
#=============================
set lmargin at screen 0.15
set rmargin at screen 0.80
set bmargin at screen 0.11
set tmargin at screen 0.21

set title ""
set xlabel "$xlabel"  font "Arial, 14"
set xtics
unset x2tics
set ylabel "$ylabel_plot4" font "Arial, 12" rotate by 0 offset 2.5,1.5
set xtics font 'Times-Roman, 14'
set ytics font 'Times-Roman, 14'
set grid noxtics ytics
set style data linespoints
set autoscale x
set autoscale y
set yrange[0:]
set ytics $ytic4
set xrange[-1:$N_record]
unset key

plot \
    $plotOption4
unset multiplot

EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Histogram image output to $pdffile"
#             if [ -s "$pdffile" ]; then
#                 rm -f "$outfile"
#             fi
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
xlabel="Pfam family (sorted in descending order of the number of homologous pairs)"
ylabel="Relative frequency (%)"
ylabel_plot2="Number of\npairs"
ylabel_plot3="Number of\npredicted\nTM proteins"
ylabel_plot4="Total\nnumber of\nproteins"
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
            -color|--color) isColor=1;shift;;
            -title|--title) title="$2";shift;;
            -multiplot|--multiplot) isMultiplot=1;;
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


stemname=${dataFile%.*}
dataFile1=$stemname.1.txt

CreateDataFile1 "$dataFile" "$dataFile1"


for isTwoClass in 0 1; do
    for isNormalize in 0 1 ; do
        Makeplot_mtp $dataFile1 $isNormalize $isTwoClass
    done
done

