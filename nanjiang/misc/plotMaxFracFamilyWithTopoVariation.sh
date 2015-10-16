#!/bin/bash
# Plot the fraction of families with topology variation cutting at different
# theshold when sorted by 
# e.g. number of pairs selected for each family, number of proteins selected
# for each family, 
# #topology comparison of four states
# #

usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle png|term|eps  set the output, default = terminal
    -outpath DIR            set the output path, (default: the same as the input file)
    -h|--help               print this help message and exit

Created 2013-11-25, updated 2013-11-25, Nanjiang Shu 
"
PlotSingle(){ #{{{
    local dataFile="$1"
    local outputSetting=
    local basename=`basename "$dataFile"`
    local outfile=
    local pdffile=

    local tag=${dataFile%.*}
    tag=${tag##*.}

    if [[ "$dataFile" =~ "_win" ]]; then
        local xlabel="${tag#sortby_}"
    else
        local xlabel="Minimal ${tag#sortby_}"
    fi

    outputSetting=
    case $outputStyle in
        term*)
            outputSetting=
            ;;
        png)
            outfile=$outpath/$basename.png
            outputSetting="
            set terminal png enhanced
            set output '$outfile' 
            "
            ;;
        eps)
            outfile=$outpath/$basename.eps
            pdffile=$outpath/$basename.pdf
            outputSetting="
            set term postscript eps color
            set output '$outfile' 
            "
            ;;
    esac

/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 pt 1 ps 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set tics font "Arial, 16"
set key top center outside
#set style data lines
#set style data points
set logscale x
set yrange [0:100]
set grid y

plot '$dataFile' using 2:3 ti '$tag' ls 1
EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Histogram image output to $pdffile"
            ;;
        *) echo "Histogram image output to $outfile" ;;
    esac

}
#}}}
PlotSingle_each_class(){ #{{{
    local dataFile="$1"
    local outputSetting=
    local basename=`basename "$dataFile"`
    local outfile=
    local pdffile=

    local tag=${dataFile%.*}
    tag=${tag##*.}
    basename=$basename.class

    if [[ "$dataFile" =~ "_win" ]]; then
        local xlabel="${tag#sortby_}"
    else
        local xlabel="Minimal ${tag#sortby_}"
    fi

    outputSetting=
    case $outputStyle in
        term*)
            outputSetting=
            ;;
        png)
            outfile=$outpath/$basename.png
            outputSetting="
            set terminal png enhanced
            set output '$outfile' 
            "
            ;;
        eps)
            outfile=$outpath/$basename.eps
            pdffile=$outpath/$basename.pdf
            outputSetting="
            set term postscript eps color
            set output '$outfile' 
            "
            ;;
    esac

/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 pt 1 ps 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

# color lines
set style line 21  lt 1 lw 1 lc rgb "red"    # IDT
set style line 22  lt 1 lw 1 lc rgb "blue"   # INV
set style line 23  lt 1 lw 1 lc rgb "yellow" # DUP
set style line 24  lt 1 lw 1 lc rgb "green"  # TM2GAP
set style line 25  lt 1 lw 1 lc rgb "violet" # Mixed
set style line 26  lt 1 lw 1 lc rgb "cyan"   # TM2SEQ
set style line 27  lt 1 lw 1 lc rgb "orange" # TM2SP

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set tics font "Arial, 16"
set key top center outside
#set style data lines
#set style data points
set logscale x
set yrange [0:100]
set grid y

plot '$dataFile' using 2:4 ti 'Inverted, $tag'    ls 22, \
    ''           using 2:5 ti 'Duplicated, $tag'  ls 23, \
    ''           using 2:6 ti 'TM2GAP, $tag'      ls 24, \
    ''           using 2:9 ti 'Mixed, $tag'       ls 25, \
    ''           using 2:7 ti 'TM2SEQ, $tag'      ls 26, \
    ''           using 2:8 ti 'TM2SP, $tag'       ls 27
EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Histogram image output to $pdffile"
            ;;
        *) echo "Histogram image output to $outfile" ;;
    esac

}
#}}}
PlotMultiple(){ #{{{
    local outputSetting=
    local basename=`basename "${dataFileList[0]}"`
    basename=$basename.mtp
    local outfile=
    local pdffile=


    local xlabel="Top N"

    outputSetting=
    case $outputStyle in
        term*)
            outputSetting=
            ;;
        png)
            outfile=$outpath/$basename.png
            outputSetting="
            set terminal png enhanced
            set output '$outfile' 
            "
            ;;
        eps)
            outfile=$outpath/$basename.eps
            pdffile=$outpath/$basename.pdf
            outputSetting="
            set term postscript eps color
            set output '$outfile' 
            "
            ;;
    esac
    local plotOption=
    local file=
    local tag=
    local cls=
    local j=
    for ((i=0;i<numFile;i++));do
        ((cls=i+1))
        ((j=i+1))
        file=${dataFileList[$i]}
        tag=${file%.*}
        tag=${tag##*.}
        plotOption="$plotOption '$file' using 1:3 ti \"$tag\" ls $cls"
        if [ $j -lt $numFile ];then
            plotOption="$plotOption,"
        fi
    done
    echo "mtp:"
    echo "$plotOption"

/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "blue" 
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set tics font "Arial, 16"
set key top center outside
set style data lines
set yrange [0:100]
set grid y

plot $plotOption
EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Histogram image output to $pdffile"
            ;;
        *) echo "Histogram image output to $outfile" ;;
    esac

}
#}}}
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
dataFileList=()
isNonOptionArg=false
xlabel="Top N"
ylabel="Frequency of families with topology variations (%)"

while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        dataFileList+=("$1")
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -outstyle|--outstyle) outputStyle=$2;shift;;
            -xlabel|--xlabel) xlabel="$2";shift;;
            -ylabel|--ylabel) ylabel="$2";shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFileList+=("$1")
    fi
    shift
done

numFile=${#dataFileList[@]}

if [ $numFile -lt 1 ];then
    echo "NumFile < 1" >&2
    exit
fi

if [ "$outpath" == "" ];then
    outpath=`dirname "${dataFileList[0]}"`
elif [ ! -d "$outpath" ];then
    mkdir -p "$outpath"
fi

for ((i=0;i<numFile;i++));do
    file=${dataFileList[$i]}
    if [ -s "$file" ];then
        PlotSingle "$file"
        numField=`awk '{print NF}' "$file" | sort -rg | head -n 1`
        if [ "$numField" != "" ];then
            if [ $numField -gt 3 ];then
                PlotSingle_each_class "$file"
            fi
        fi
    fi
done

if [ $numFile -gt 1 ];then
    PlotMultiple
fi

