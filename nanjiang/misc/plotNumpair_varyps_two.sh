#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #


usage="
Usage: $0 datafile1 datafile2
Description: 
Options:
    -outstyle|-o png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2012-09-21, updated 2012-09-21, Nanjiang Shu 
"
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outpath=
datafileList=()
isNonOptionArg=false
xlabel="Reliability score threshold"
ylabel="Number of pairs"

while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        datafileList+=("$1")
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
        datafileList+=("$1")
    fi
    shift
done

numFile=${#datafileList[@]}
if [ $numFile -lt 2 ]; then
    echo "too few datafile (<2)"
    exit 1
fi

datafile1=${datafileList[0]}
datafile2=${datafileList[1]}



if [ ! -s "$datafile1" -o ! -s "$datafile2" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    exit
fi


if [ "$outpath" == "" ]; then
    outpath=`dirname $datafile1`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

col=`head -n 1 $datafile1 | awk '{print NF}'`

basename=`basename "$datafile1"` 
basename=$basename.two

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
set style line 1 lt 1 lw 1 lc rgb "green"
set style line 2 lt 1 lw 1 lc rgb "violet"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set key top center outside

set style data histogram 
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'

plot \
    '$datafile1' using $col:xtic(2)   title 'Topcons' ls 1 , \
    '$datafile2' using $col:xtic(2)   title 'Topcons_single' ls 2
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


