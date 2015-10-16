#!/bin/bash
# Filename:  plotHelixNumContUnmappedTM_cmpdup.sh
# Description:
#   Draw histogram for the number of TM helices of indels

# Groups for pairwise_comparison_method 3:
#   0. TM2GAP in all
#   1. TM2SEQ in all
#   2. TM2GAP in TM2GAP
#   3. TM2GAP in Mixed
#   4. TM2GAP in DUP
#   5. TM2SEQ in TM2SEQ
#   6. TM2SEQ in Mixed

usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle,-o STR    Set the output, (default: eps)
                        can be png, term, eps
    -outpath DIR        Set the output path, (default: ./)
    -logscale y|n       Plot y axis as logscale, (default: yes)
    -xlabel STR         Set xlabel, (default: Sequence identity)
    -ylabel STR         Set Ylabel
    -h,--help           Print this help message and exit

Created 2013-10-28, updated 2013-10-28, Nanjiang Shu
"
cmd="$*"
PrintHelp(){
    echo "$usage"
}
MakePlot() { #{{{  # using linespoints
    local mode=$1
    local b1=$basename
    local plotOption=""
    case $mode in 
        0|1|2|3|4|5)
            li=(1  9)
            d1=${li[0]}
            ((addon=mode*8))
            ((col1=addon+1))
            ((col2=addon+2))
            plotOption="'$dataFile' using $col1:$col2 ti 'pos 0' ls 1"
            num=${#li[@]}
            for((i=0;i<num;i++));do
                d1=${li[$i]}
                ((col1=addon+d1*2+1))
                ((col2=addon+d1*2+2))
                x1=0.$d1
                ((j=i+1))
                if [ $j -lt $num ]; then 
                    d2=${li[$j]}
                    x2=0.$d2
                else
                    x2=1.0
                fi
                ((ls=i+2))
                plotOption="$plotOption, '' using $col1:$col2 ti 'pos $x1 - $x2' ls $ls"
            done
            b1=$b1.$mode
            ;;
        45-all)
            plotOption="'$dataFile' using 109:110 ti 'TM2GAP' ls 1,\
                '' using 131:132 ti 'TM2SEQ' ls 2"
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
    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set xrange[1:14]
#set key top center outside
set datafile missing "?"
set style data linespoints 
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
MakePlot1() { #{{{  #using points
    local mode=$1
    local b1=$basename
    local plotOption=""
    case $mode in 
        0|1|2|3|4|5)
            d1=${li[0]}
            ((addon=mode*8))
            ((col1=addon+1))
            ((col2=addon+2))
            plotOption="'$dataFile' using $col1:$col2 ti 'N-terminal' ls 1"
            ((col1=addon+3))
            ((col2=addon+4))
            plotOption="$plotOption, '' using $col1:$col2 ti 'Internal' ls 2"
            ((col1=addon+5))
            ((col2=addon+6))
            plotOption="$plotOption, '' using $col1:$col2 ti 'C-terminal' ls 3"
            b1=$b1.$mode
            ;;
        45-all)
            m=4
            ((addon=m*8))
            ((col1=addon+7))
            ((col2=addon+8))
            plotOption="'$dataFile' using $col1:$col2 ti 'TM2GAP' ls 1"
            m=5
            ((addon=m*8))
            ((col1=addon+7))
            ((col2=addon+8))
            plotOption="$plotOption, '' using $col1:$col2 ti 'TM2SEQ' ls 2"
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
    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 3 lw 1 pt 5 ps 2 lc rgb "red"
set style line 2 lt 3 lw 1 pt 6 ps 2 lc rgb "green"
set style line 3 lt 3 lw 1 pt 7 ps 2 lc rgb "blue" 
set style line 4 lt 3 lw 1 pt 8 ps 2 lc rgb "violet"
set style line 5 lt 3 lw 1 pt 9 ps 2 lc rgb "cyan" 
set style line 6 lt 3 lw 1 pt 10 ps 2 lc rgb "orange"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set xrange[1:14]
#set key top center outside
set datafile missing "?"
#set style data linespoints 
set style data points 
set style fill solid
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
MakePlot2() { #{{{ # using lines points, dashed line, black and white, 
    local mode=$1
    local b1=$basename
    local plotOption=
    local addon=
    local col1=
    local col2=
    case $mode in 
        0|1|2|3|4|5|6)
            d1=${li[0]}
            ((addon=mode*10))
            ((col1=addon+1))
            ((col2=addon+2))
            plotOption="'$dataFile' using $col2:xtic($col1) ti 'N-terminus' ls 1"
            ((col1=addon+3))
            ((col2=addon+4))
            plotOption="$plotOption, '' using $col2:xtic($col1) ti 'Internal' ls 2"
            ((col1=addon+5))
            ((col2=addon+6))
            plotOption="$plotOption, '' using $col2:xtic($col1) ti 'C-terminus' ls 3"
            ((col1=addon+7))
            ((col2=addon+8))
            plotOption="$plotOption, '' using $col2:xtic($col1) ti 'Duplicated' ls 4"
            b1=$b1.$mode
            ;;
        0-all)
            m=0
            ((addon=m*10))
            ((col1=addon+9))
            ((col2=addon+10))
            plotOption="'$dataFile' using $col2:xtic($col1) ti 'TM2GAP in all' ls 0"
            b1=$b1.$mode
            ;;
        1-all)
            m=1
            ((addon=m*10))
            ((col1=addon+9))
            ((col2=addon+10))
            plotOption="'$dataFile' using $col2:xtic($col1) ti 'TM2SEQ in all' ls 1"
            b1=$b1.$mode
            ;;
        4-all)
            m=4
            ((addon=m*10))
            ((col1=addon+9))
            ((col2=addon+10))
            plotOption="'$dataFile' using $col2:xtic($col1) ti 'TM2GAP in DUP' ls 4"
            b1=$b1.$mode
            ;;
        5-all)
            m=5
            ((addon=m*10))
            ((col1=addon+9))
            ((col2=addon+10))
            plotOption="'$dataFile' using $col2:xtic($col1) ti 'TM2SEQ in TM2SEQ' ls 5"
            b1=$b1.$mode
            ;;
        01-all)
            m=0
            ((addon=m*10))
            ((col1=addon+9))
            ((col2=addon+10))
            plotOption="'$dataFile' using $col2:xtic($col1) ti 'TM2GAP in all' ls 0"
            m=1
            ((addon=m*10))
            ((col1=addon+9))
            ((col2=addon+10))
            plotOption="$plotOption, '' using $col2:xtic($col1) ti 'TM2SEQ in all' ls 1"
            b1=$b1.$mode
            ;;
        *)return 1
            ;;
    esac

    #echo "$plotOption"
    local ylabel=$ylabel
    local logscaleopt=
    if [ $isLogScale -eq 1 ];then
        b1=$b1.logscale
        ylabel="$ylabel (log scale)"
        logscaleopt="set logscale y"
    fi

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
set style line 1 lt 3 lw 1 pt 7 ps 2 lc rgb "black"
set style line 2 lt 3 lw 1 pt 8 ps 2 lc rgb "black"
set style line 3 lt 3 lw 1 pt 9 ps 2 lc rgb "black" 
set style line 4 lt 2 lw 3 pt 3 ps 2 lc rgb "black"
set style line 5 lt 3 lw 1 pt 8 ps 2 lc rgb "black" 

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set xrange[0:13]
$logscaleopt
#set key top center outside
set datafile missing "?"
set style data linespoints
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
xlabel="Number of TM helices"
ylabel="Occurrences"
isNonOptionArg=false
isLogScale=1

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
            -logscale|--logscale)
                case ${2:0:1} in 
                    y|Y) isLogScale=1;;
                    *) isLogScale=0;;
                esac
                shift;;
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

# MakePlot "everything"
# MakePlot "all"
# MakePlot "both"
# MakePlot "only"
for item in "0" "1" "2" "3" "4" "5" "0-all" "1-all" "4-all" "5-all" "45-all" "01-all"; do
    MakePlot2 $item
done
