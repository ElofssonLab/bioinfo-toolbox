#!/bin/bash

usage="
Usage: $0 -tm2tm FILE -tm2gap FILE -tm2seq FILE -tm2sp FILE
Description: 
Options:
    -o, -outstyle STR   set the output, (default: eps)
                        can be png, term, eps
    -minseqidt  FLOAT
    -maxseqidt  FLOAT
    -outname      STR   set the name of the outfile
    -h, --help          print this help message and exit

Created 2015-04-20, updated 2015-04-20, Nanjiang Shu 

Example:
"
PrintHelp(){
    echo "$usage"
}

Makeplot(){ #{{{  #not show number of cases
    local outputSetting=""
    local outfile=
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$outname.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$outname.eps
        pdffile=$outname.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac
    pdffile=$outname.pdf
    plotOption="\
        '$tm2tm' using 12:13 ti 'TM2TM' ls 1, \
        '$tm2gap' using 12:(\$13+10) ti 'TM2GAP' ls 2,\
        '$tm2seq' using 12:13 ti 'TM2SEQ' ls 3, \
        '$tm2sp' using 12:13 ti 'TM2SP' ls 4
        "
    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"
set title "$title" font "Arial, 20"
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"
set xzeroaxis
set yzeroaxis

plot \
    $plotOption
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

title=
outputStyle=eps
outpath=
xlabel="DeltaG"
ylabel="DeltaG"
ylabel_plot2="Number of pairs"
isNonOptionArg=false
minSEQIDT=0
maxSEQIDT=100

tm2sp=
tm2tm=
tm2gap=
tm2seq=

fs1="fs solid"
fs2="fs pattern 2"
fs3="fs pattern 3"
fs4="fs pattern 1"
fs5="fs pattern 4"
fs6="fs pattern 7"

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
            -title|--title) title="$2";shift;;
            -tm2tm|--tm2tm) tm2tm="$2";shift;;
            -tm2gap|--tm2gap) tm2gap="$2";shift;;
            -tm2seq|--tm2seq) tm2seq="$2";shift;;
            -tm2sp|--tm2sp) tm2sp="$2";shift;;
            -minseqidt)minSEQIDT=$2;shift;;
            -maxseqidt)maxSEQIDT=$2;shift;;
            -outname|--outname) outname=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFile=$1
    fi
    shift
done

if [ "$outname" == "" ];then
    echo "outname not set exit"
    exit 1
fi
outpath=`dirname $outname`
if [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

tmpdir=$(mktemp -d /tmp/tmpdir.plotTMHelixSegmentDGScoreXY.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }

if [ "$minSEQIDT" != "" -o "$maxSEQIDT" != "" ];then
    awk -v v1=$minSEQIDT -v v2=$maxSEQIDT '{if($4>=v1&&$4<v2)print }' $tm2tm > $tmpdir/tm2tm.txt
    awk -v v1=$minSEQIDT -v v2=$maxSEQIDT '{if($4>=v1&&$4<v2)print }' $tm2gap > $tmpdir/tm2gap.txt
    awk -v v1=$minSEQIDT -v v2=$maxSEQIDT '{if($4>=v1&&$4<v2)print }' $tm2sp > $tmpdir/tm2sp.txt
    awk -v v1=$minSEQIDT -v v2=$maxSEQIDT '{if($4>=v1&&$4<v2)print }' $tm2seq > $tmpdir/tm2seq.txt
    tm2tm=$tmpdir/tm2tm.txt
    tm2sp=$tmpdir/tm2sp.txt
    tm2seq=$tmpdir/tm2seq.txt
    tm2gap=$tmpdir/tm2gap.txt
    outname=$outname.$minSEQIDT-$maxSEQIDT
    title="seqidt $minSEQIDT - $maxSEQIDT"
fi

Makeplot

rm -rf $tmpdir
