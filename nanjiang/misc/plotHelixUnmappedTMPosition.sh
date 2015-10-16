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

Created 2012-09-25, updated 2012-11-19, Nanjiang Shu 
"
cmd="$*"
PrintHelp(){
    echo "$usage"
}
MakePlot(){ #{{{ # using lines points
    local mode=$1
    local b1=$basename
    local plotOption=""
    if [ "$mode" == "everything" ]; then
        plotOption="\
        '$dataFile' using 1:2  title 'TM2GAP, Only TM2GAP' ls 1,\
        ''          using 3:4  title 'TM2SEQ, Both' ls 2,\
        ''          using 5:6  title 'TM2SEQ, Only TM2SEQ' ls 3,\
        ''          using 7:8  title 'TM2SEQ, Both' ls 4,\
        ''          using 9:10  title 'TM2SEQ, All' ls 5,\
        ''          using 11:12  title 'TM2SEQ, All' ls 6
        "
        b1=$b1.everything
    elif [ "$mode" == "all" ]; then
        plotOption="\
        '$dataFile' using 9:10  title 'TM2GAP, All' ls 5,\
        ''          using 11:12  title 'TM2SEQ, All' ls 6
        "
        b1=$b1.all
    elif [ "$mode" == "both" ]; then
        plotOption="\
        '$dataFile' using 3:4  title 'TM2GAP, Both' ls 2,\
        ''          using 7:8  title 'TM2SEQ, Both' ls 4
        "
        b1=$b1.both
    elif [ "$mode" == "only" ]; then
        plotOption="\
        '$dataFile' using 1:2  title 'TM2GAP, Only TM2GAP' ls 1,\
        ''          using 5:6  title 'TM2SEQ, Only TM2SEQ' ls 3
        "
        b1=$b1.only
    fi
    
    #echo "$plotOption"

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
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
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
MakePlot1(){ #{{{ # using histogram
    local mode=$1
    local b1=$basename
    local plotOption=""
    local tmpdatafile=$(mktemp /tmp/tmp.plotHelixUnmappedTMPosition.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    trap 'rm -f "$tmpdatafile"' INT TERM EXIT
    msg=`grep N-terminal $dataFile`
    if [  "$msg" == "" ]; then
        cat $dataFile | sed 's/0*[ \t]/ /g' | sed 's/\.[ \t]/ /g' > $tmpdatafile
    else
        cat $dataFile  > $tmpdatafile
    fi


    case $mode in 
        everything)
            plotOption="\
                '$tmpdatafile' using 2:xtic(1) title 'TM2GAP, Only' fs pattern 1 ls 1,\
                ''          using 4:xtic(3) title 'TM2GAP, Both' fs pattern 1 ls 2,\
                ''          using 6:xtic(5)  title 'TM2SEQ, Only' fs patten 2 ls 3,\
                ''          using 8:xtic(7)  title 'TM2SEQ, Both' fs pattern 2 ls 4,\
                ''          using 10:xtic(9)  title 'TM2GAP, All' fs solid 1 ls 5,\
                ''          using 12:xtic(11)  title 'TM2SEQ, All' fs solid 1 ls 6
            "
            b1=$b1.everything
            ;;
        all)
            plotOption="\
                '$tmpdatafile' using 10:xtic(9)  title 'TM2GAP' fs solid 1 ls 5,\
                ''          using 12:xtic(11)  title 'TM2SEQ' fs pattern 1 ls 6
            "
            b1=$b1.all
            ;;
        both)
            plotOption="\
                '$tmpdatafile' using 4:xtic(3)  title 'TM2GAP' fs solid 1 ls 2,\
                ''          using 8:xtic(7)  title 'TM2SEQ' fs pattern 1 ls 4
            "
            b1=$b1.both
            ;;
        only)
            plotOption="\
                '$tmpdatafile' using 2:xtic(1)  title 'TM2GAP' fs solid 1 ls 1,\
                ''          using 6:xtic(5)  title 'TM2SEQ' fs pattern 1 ls 3
            "
            b1=$b1.only
            ;;
        0)
            plotOption="\
                '$tmpdatafile' using 2:xtic(1)  title 'TM2GAP, Only' fs solid 1 ls 1
            "
            b1=$b1.$mode
            ;;
        1)
            plotOption="\
                '$tmpdatafile' using 4:xtic(3)  title 'TM2GAP, Both' fs solid 1 ls 1
            "
            b1=$b1.$mode
            ;;
        2)
            plotOption="\
                '$tmpdatafile' using 6:xtic(5)  title 'TM2SEQ, Only' fs solid 1 ls 1
            "
            b1=$b1.$mode
            ;;
        3)
            plotOption="\
                '$tmpdatafile' using 8:xtic(7)  title 'TM2SEQ, Both' fs solid 1 ls 1
            "
            b1=$b1.$mode
            ;;
        4)
            plotOption="\
                '$tmpdatafile' using 10:xtic(9)  title 'TM2GAP, All' fs solid 1 ls 1
            "
            b1=$b1.$mode
            ;;
        5)
            plotOption="\
                '$tmpdatafile' using 12:xtic(11)  title 'TM2SEQ, All' fs solid 1 ls 1
            "
            b1=$b1.$mode
            ;;
    esac
    
    #echo "$plotOption"

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
set style line 1 lt 1 lw 1 lc rgb "black"
set style line 2 lt 1 lw 1 lc rgb "grey70"
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
xlabel="Normalized sequence position"
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
case $dataFile in 
    *3bin*)
        xlabel="Sequence position"
        ;;
    *)
        xlabel="Normalized sequence position"
        ;;
esac

#MakePlot1 "everything"
for item in all both only 0 1 2 3 4 5; do 
    MakePlot1 $item
done

