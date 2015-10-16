#!/bin/bash
# Filename: plotHelixNumContUnmappedTM.withHistogram.sh 
# Description:
#   Draw histogram multiplot

# multiplot
usage="
Usage: $0 numcontHistFile  unmappedTMpositionHistFile
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
MakePlot_multiplot() { #{{{ # using lines points, dashed line, black and white, 
    local mode=$1
    local b1=$basename
    local plotOption1=
    local plotOption2=
    local addon=
    local col1=
    local col2=
    local col3=
    local title2=
    case $mode in 
        0|1|2|3|4|5|6)
            d1=${li[0]}
            ((addon=mode*10))
            ((col1=addon+1))
            ((col2=addon+2))
            plotOption1="'$dataFile1' using $col2:xtic($col1) ti 'N-terminus' ls 1"
            ((col1=addon+3))
            ((col2=addon+4))
            plotOption1="$plotOption1, '' using $col2:xtic($col1) ti 'Internal' ls 2"
            ((col1=addon+5))
            ((col2=addon+6))
            plotOption1="$plotOption1, '' using $col2:xtic($col1) ti 'C-terminus' ls 3"
            ((col1=addon+7))
            ((col2=addon+8))
            plotOption1="$plotOption1, '' using $col2:xtic($col1) ti 'Duplicated' ls 4"

            ((addon=mode*3))
            ((col1=addon+2))
            ((col2=addon+3))
            ((col3=addon+4))
            plotOption2="\
                '$dataFile2' using $col1:xtic(1)  title '1 TM' fs solid 1 ls 11, \
                ''           using $col2:xtic(1)  title '2 TM' fs solid 1 ls 12, \
                ''           using $col3:xtic(1)  title '>2 TM' fs solid 1 ls 13 \
                "

            if [ "$mode" == "0" ];then
                title2="TM2GAP"
            elif [ "$mode" == "1" ];then
                title2="TM2SEQ"
            elif [ "$mode" == "2" ];then
                title2="TM2GAP in TM2GAP"
            elif [ "$mode" == "3" ];then
                title2="TM2GAP in Mixed"
            elif [ "$mode" == "4" ];then
                title2="TM2GAP in DUP"
            elif [ "$mode" == "5" ];then
                title2="TM2SEQ in TM2SEQ"
            elif [ "$mode" == "6" ];then
                title2="TM2SEQ in Mixed"
            fi
            b1=$b1.$mode.multiplot
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

set style line 11 lt 1 lw 1 pt 1 ps 1 lc rgb "grey90"
set style line 12 lt 1 lw 1 pt 1 ps 1 lc rgb "grey50"
set style line 13 lt 1 lw 1 pt 1 ps 1 lc rgb "black"

set style line 1 lt 3 lw 1 pt 7 ps 2 lc rgb "black"
set style line 2 lt 3 lw 1 pt 8 ps 2 lc rgb "black"
set style line 3 lt 3 lw 1 pt 9 ps 2 lc rgb "black" 
set style line 4 lt 2 lw 3 pt 3 ps 2 lc rgb "black"
set style line 5 lt 3 lw 1 pt 8 ps 2 lc rgb "black" 

set multiplot
set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" offset 1 font "Arial, 16"
set xrange[0:13]
$logscaleopt
#set key top center outside
set datafile missing "?"
set style data linespoints
set tics font 'Times-Roman,16'

plot $plotOption1

# the smaller plot
set title "$title3"
set size 0.45,0.45
set origin 0.5,0.4
set autoscale x
set autoscale y
unset logscale y
#set xlabel "Location"  font "Arial, 12"
unset xlabel
set ylabel "Frequency" offset 3 font "Arial, 12"
set datafile missing "?"
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set tics font 'Times-Roman,12'
set xtics rotate by -30
plot $plotOption2

unset multiplot

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
MakePlot_loc() { #{{{ # using lines points, dashed line, black and white, 
    # N-term, C-term and internal
    local mode=$1
    local b1=$basename
    b1=$b1.loc
    local plotOption1=
    local plotOption2=
    local addon=
    local col1=
    local col2=
    local col3=
    local title2=
    case $mode in 
        0|1|2|3|4|5|6)
            d1=${li[0]}
            ((addon=mode*10))
            ((col1=addon+1))
            ((col2=addon+2))
            plotOption1="'$dataFile1' using $col2:xtic($col1) ti 'N-terminus' ls 1"
            ((col1=addon+3))
            ((col2=addon+4))
            plotOption1="$plotOption1, '' using $col2:xtic($col1) ti 'Internal' ls 2"
            ((col1=addon+5))
            ((col2=addon+6))
            plotOption1="$plotOption1, '' using $col2:xtic($col1) ti 'C-terminus' ls 3"

            if [ "$mode" == "0" ];then
                title2="TM2GAP"
            elif [ "$mode" == "1" ];then
                title2="TM2SEQ"
            elif [ "$mode" == "2" ];then
                title2="TM2GAP in TM2GAP"
            elif [ "$mode" == "3" ];then
                title2="TM2GAP in Mixed"
            elif [ "$mode" == "4" ];then
                title2="TM2GAP in DUP"
            elif [ "$mode" == "5" ];then
                title2="TM2SEQ in TM2SEQ"
            elif [ "$mode" == "6" ];then
                title2="TM2SEQ in Mixed"
            fi
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

set style line 11 lt 1 lw 1 pt 1 ps 1 lc rgb "grey90"
set style line 12 lt 1 lw 1 pt 1 ps 1 lc rgb "grey50"
set style line 13 lt 1 lw 1 pt 1 ps 1 lc rgb "black"

set style line 1 lt 2 lw 2 pt 5 ps 2 lc rgb "red"
set style line 2 lt 2 lw 2 pt 7 ps 2 lc rgb "grey20"
set style line 3 lt 2 lw 2 pt 9 ps 2 lc rgb "blue" 

set style line 4 lt 2 lw 3 pt 3 ps 2 lc rgb "black"
set style line 5 lt 3 lw 1 pt 8 ps 2 lc rgb "black" 

set multiplot
set title ""
set xlabel "$xlabel"  font "$font_label"
set ylabel "$ylabel" offset -1 font "$font_label"
set xrange[0:9]
set yrange[-1:]
$logscaleopt
set key top right font "$font_key" spacing 1.5
#set datafile missing "?"

set style data linespoints
# set style data histogram
# set style histogram cluster gap 1
# set style fill solid border -1
# set boxwidth 0.9

set tics font "$font_tics"

plot $plotOption1

EOF

case $outputStyle in
    eps) my_epstopdf $outfile
        echo "Histogram image output to $pdffile"
#         if [ -s "$pdffile" ]; then
#             rm -f "$outfile"
#         fi
        ;;
    *) echo "Histogram image output to $outfile" ;;
esac

}
#}}}
MakePlot_duplication() { #{{{ # using lines points, dashed line, black and white, 
    local mode=$1
    local b1=$basename
    b1=$b1.dup
    local plotOption1=
    local plotOption2=
    local addon=
    local col1=
    local col2=
    local col3=
    local title2=
    case $mode in 
        0|1|2|3|4|5|6)
            d1=${li[0]}
            ((addon=mode*10))
            ((col1=addon+7))
            ((col2=addon+8))
            plotOption1="'$dataFile1' using $col2:xtic($col1) ti 'Duplicated' ls 4"

            if [ "$mode" == "0" ];then
                title2="TM2GAP"
            elif [ "$mode" == "1" ];then
                title2="TM2SEQ"
            elif [ "$mode" == "2" ];then
                title2="TM2GAP in TM2GAP"
            elif [ "$mode" == "3" ];then
                title2="TM2GAP in Mixed"
            elif [ "$mode" == "4" ];then
                title2="TM2GAP in DUP"
            elif [ "$mode" == "5" ];then
                title2="TM2SEQ in TM2SEQ"
            elif [ "$mode" == "6" ];then
                title2="TM2SEQ in Mixed"
            fi
            b1=$b1.$mode
            ;;
        *)return 1
            ;;
    esac
    echo "plotOption1=$plotOption1"

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

set style line 11 lt 1 lw 1 pt 1 ps 1 lc rgb "grey90"
set style line 12 lt 1 lw 1 pt 1 ps 1 lc rgb "grey50"
set style line 13 lt 1 lw 1 pt 1 ps 1 lc rgb "black"

set style line 1 lt 3 lw 1 pt 7 ps 2 lc rgb "black"
set style line 2 lt 3 lw 1 pt 8 ps 2 lc rgb "black"
set style line 3 lt 3 lw 1 pt 9 ps 2 lc rgb "black" 
set style line 4 lt 1 lw 1 pt 3 ps 2 lc rgb "grey"
set style line 5 lt 3 lw 1 pt 8 ps 2 lc rgb "black" 

set title ""
set xlabel "$xlabel"  font "$font_label"
set ylabel "$ylabel" offset -1 font "$font_label"
set key top right font "$font_key" spacing 1.5
set tics font "$font_tics"
set xrange[-1:9]
#set yrange[1:400]
$logscaleopt
set datafile missing "?"

#set style data linespoints

set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9


plot $plotOption1
EOF

case $outputStyle in
    eps) my_epstopdf $outfile
        echo "Histogram image output to $pdffile"
#         if [ -s "$pdffile" ]; then
#             rm -f "$outfile"
#         fi
        ;;
    *) echo "Histogram image output to $outfile" ;;
esac

}
#}}}
MakePlot_histogram() { #{{{ # using lines points, dashed line, black and white, 
    local mode=$1
    local b1=$basename
    local plotOption1=
    local plotOption2=
    local addon=
    local col1=
    local col2=
    local col3=
    local title2=
    b1=$b1.hist
    case $mode in 
        0|1|2|3|4|5|6)
            d1=${li[0]}
            ((addon=mode*3))
            ((col1=addon+2))
            ((col2=addon+3))
            ((col3=addon+4))
            plotOption2="\
                '$dataFile2' using $col1:xtic(1)  title '1 TM' fs solid 1 ls 11, \
                ''           using $col2:xtic(1)  title '2 TM' fs solid 1 ls 12, \
                ''           using $col3:xtic(1)  title '>2 TM' fs solid 1 ls 13 \
                "

            if [ "$mode" == "0" ];then
                title2="TM2GAP"
            elif [ "$mode" == "1" ];then
                title2="TM2SEQ"
            elif [ "$mode" == "2" ];then
                title2="TM2GAP in TM2GAP"
            elif [ "$mode" == "3" ];then
                title2="TM2GAP in Mixed"
            elif [ "$mode" == "4" ];then
                title2="TM2GAP in DUP"
            elif [ "$mode" == "5" ];then
                title2="TM2SEQ in TM2SEQ"
            elif [ "$mode" == "6" ];then
                title2="TM2SEQ in Mixed"
            fi
            b1=$b1.$mode
            ;;
        *)return 1
            ;;
    esac

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

set style line 11 lt 1 lw 1 pt 1 ps 1 lc rgb "grey90"
set style line 12 lt 1 lw 1 pt 1 ps 1 lc rgb "grey50"
set style line 13 lt 1 lw 1 pt 1 ps 1 lc rgb "black"

set style line 1 lt 3 lw 1 pt 7 ps 2 lc rgb "black"
set style line 2 lt 3 lw 1 pt 8 ps 2 lc rgb "black"
set style line 3 lt 3 lw 1 pt 9 ps 2 lc rgb "black" 
set style line 4 lt 2 lw 3 pt 3 ps 2 lc rgb "black"
set style line 5 lt 3 lw 1 pt 8 ps 2 lc rgb "black" 


set title ""
set autoscale x
set autoscale y
unset logscale y
unset xlabel
set ylabel "Frequency" offset -1 font "$font_label"
set key top right font "$font_key" spacing 1.5
set datafile missing "?"
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set tics font "$font_tics"
#set xtics rotate by -30
plot $plotOption2
EOF

case $outputStyle in
    eps) my_epstopdf $outfile
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

outputStyle=eps
outpath=
dataFileList=()
xlabel="Number of TM helices"
ylabel="Occurrences"
isNonOptionArg=false
isLogScale=1

while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        dataFileList+=("$1")
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
        dataFileList+=("$1")
    fi
    shift
done


numFile=${#dataFileList[@]}
if [ $numFile -ne 2 ];then
    echo "number of input files != 2. exit" >&2
    exit
fi

dataFile1=${dataFileList[0]}
dataFile2=${dataFileList[1]}

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile1`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

if [ ! -f "$dataFile1" ]; then 
    echo "Error! dataFile1 = \"$dataFile1\" does not exist. Exit..."
    echo "$cmd"
    exit
fi
if [ ! -f "$dataFile2" ]; then 
    echo "Error! dataFile2 = \"$dataFile2\" does not exist. Exit..."
    echo "$cmd"
    exit
fi

font_label="Arial, 24"
font_tics="Arial, 24"
font_key="Arial, 24"

basename=`basename "$dataFile2"` 
dirname=`dirname $dataFile2`

# MakePlot "everything"
# MakePlot "all"
# MakePlot "both"
# MakePlot "only"
for item in "0" "1" "2" "3" "4" "5" "6"; do
#for item in 0; do
    MakePlot_multiplot $item
    MakePlot_histogram $item
    MakePlot_duplication $item
    MakePlot_loc $item
done
