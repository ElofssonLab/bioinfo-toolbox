#!/bin/bash
# draw the histogram given the text file output by calc_dG.pl
# #
progname=`dirname $0`
usage="
Usage: $progname datafile

Description:
    Draw the histogram given the text file output by calc_dG.pl

Options:
    -outstyle|-o png|term|eps : set the output, default = terminal
    -outpath <path>           : set the output path, default =./
    -h|--help                 : print this help message and exit

Created 2010-09-03, updated 2010-09-03, Nanjiang
"
PrintHelp() {
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=term
outpath=
dataFile=
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

basename=`basename "$dataFile"` 

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
set term postscript eps enhanced
set output '$outfile' 
    "
    ;;
esac

tmpfile1=`tempfile`
awk '{if (NR>2)print $2}'  $dataFile  > $tmpfile1
mean="`stat.pl -i $tmpfile1 | awk '{print $4}'` "
rm -f $tmpfile1

tmpfile=`tempfile`
trap 'rm -f "$tmpfile"' INT TERM EXIT
awk '{if (NR>2)print $2}'  $dataFile | histogram > $tmpfile

/usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt rgb "black" lw 1    #DIFF
set style line 2 lt rgb "cyan" lw 1     #SHIFT
set style line 3 lt rgb "violet" lw 1     #INV
set style line 4 lt rgb "green" lw 1    #OK

set title "DG histogram for TM fragments"
#set label "$st" center upper
set xlabel "DG values (mean = $mean)" #0.0,-0.8 # move the xlabel 1 letter downwards
set ylabel "Normalized frequence"
set style data linespoints
# set key right outside box
# set style histogram rowstacked
# set style fill solid border -1
# set boxwidth 0.8
# set xtic rotate by -45
# #set tics font "Times-Roman,8"
# #set xtics font  "Helvetica,8"
# #set xtics font 'arial,7'
# set tics font 'Verdana,7'
# set bmargin 5
plot "$tmpfile" using 1:2 title "$dataFile" 

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

rm -f $tmpfile
