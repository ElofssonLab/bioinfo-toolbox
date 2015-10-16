#!/bin/bash

usage="
Usage:    plot_cmpclass_all.sh  all.generalinfo.txt
Options:
    -outstyle  png|term|eps     set the output, (default: eps)
    -outpath   DIR              set output dir
    -ordertype STR              set order type
                                numseq:
                                idt:
                                shift:
                                inv:
                                inv_shift:
                                diff:
    -mergeOKSHIFT y|n           Merge OK to SHIFT, (default: yes)
    -min-numseq  INT            Minimal number of sequence in family
                                (default: 20)
    -min-peridt  float          Minimal percentage of identical group
                                (default: 30)
                                
    -h|--help                   print this help message and exit

Created 2012-04-03, updated 2012-04-03, Nanjiang Shu 
Examples:
    plot_cmpclass_all.sh all.generalinfo.txt -ordertype IDT
    
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
isNonOptionArg=false
infile=all.generalinfo.txt
isMergeOK_to_SHIFT=1
ordertype=IDT
min_numseq=20
min_peridt=30.0

#NUMSEQ
#SHIFT
#INV
#INV_SHIFT
#DIFF

while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        infile=$1
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -outstyle|--outstyle) outputStyle=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -ordertype|--ordertype) ordertype=$2;shift;;
            -min-numseq|--min-numseq ) min_numseq=$2;shift;;    
            -min-peridt|--min-peridt ) min_peridt=$2;shift;;    
            -mergeOKSHIFT|--mergeOKSHIFT)
                opt=$2
                if [ "${2:0:1}" == "y" -o "${2:0:1}" == "Y" ]; then
                    isMergeOK_to_SHIFT=1
                else
                    isMergeOK_to_SHIFT=0
                fi
                shift;
                ;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

if [ ! -f "$infile" ]; then 
    echo infile $infile not exist. exit. >&2
    exit 1
fi

if [ "$outpath" == "" ]; then
    outpath=`dirname $infile`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi

basename=`basename "$infile"` 
rootname=${basename%.*}


tmpdir=$(mktemp -d /tmp/tmpdir.plot_cmpclass_all.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }  
trap 'rm -rf "$tmpdir"' INT TERM EXIT

data1=$tmpdir/data1.txt
data2=$tmpdir/data2.txt
plotOption1=
if [ $isMergeOK_to_SHIFT -eq 0 ]; then
    # numseq IDT% OK% SHIFT% INV% INV_SHIFT% DIFF%
    awk -v MIN_NUMSEQ=$min_numseq -v MIN_PERIDT=$min_peridt '
    {if($3>=MIN_NUMSEQ && $12>=MIN_PERIDT) print  $3, $12, $13, $14, $15, $16, $17}
    ' $infile > $data1
    if [ ! -s $data1 ]; then
        echo data1 $data1 is empty. abort. >&2
    else
        numPfam=`cat $data1 | wc -l`
        case $ordertype in 
            IDT)
                sort -k2,2rg -k3,3rg -k4,4rg -k5,5rg -k6,6rg -k7,7rg $data1 |\
                    awk '{print NR-1, $0}' > $data2
                plotOption1="
plot '$data2' using 8:xtic(1) ls 6 t 'DIFF', \
    ''        using 7:xtic(1) ls 5 t 'INV_SHIFT', \
     ''       using 6:xtic(1) ls 4 t 'INV', \
     ''       using 5:xtic(1) ls 3 t 'SHIFT', \
     ''       using 4:xtic(1) ls 2 t 'OK', \
     ''       using 3:xtic(1) ls 1 t 'IDT'
"
                ;;
        esac
    fi
else 
    # numseq IDT% OK+SHIFT% INV% INV_SHIFT% DIFF%
    awk -v MIN_NUMSEQ=$min_numseq -v MIN_PERIDT=$min_peridt '
    {if($3>=MIN_NUMSEQ && $12 >= MIN_PERIDT) print  $3, $12, $13+$14, $15, $16, $17}
    ' $infile > $data1
    if [ ! -s $data1 ]; then
        echo data1 $data1 is empty. abort. >&2
    else
        numPfam=`cat $data1 | wc -l`
        case $ordertype in 
            NUMSEQ|numseq)
                sort -k1,1g  $data1 |\
                    awk '{print NR-1, $0}' > $data2
                plotOption1="
plot '$data2' using 7:xtic(1) ls 6 t 'DIFF', \
    ''        using 6:xtic(1) ls 5 t 'INV_SHIFT', \
     ''       using 5:xtic(1) ls 4 t 'INV', \
     ''       using 4:xtic(1) ls 3 t 'SHIFT', \
     ''       using 3:xtic(1) ls 1 t 'IDT'
"
                ;;
            IDT|idt)
                sort -k2,2rg -k3,3rg -k4,4rg -k5,5rg -k6,6rg $data1 |\
                    awk '{print NR-1, $0}' > $data2
                plotOption1="
plot '$data2' using 7:xtic(1) ls 6 t 'DIFF', \
    ''        using 6:xtic(1) ls 5 t 'INV_SHIFT', \
     ''       using 5:xtic(1) ls 4 t 'INV', \
     ''       using 4:xtic(1) ls 3 t 'SHIFT', \
     ''       using 3:xtic(1) ls 1 t 'IDT'
"
                ;;
            SHIFT|shift)
                sort -k3,3rg -k4,4rg -k5,5rg -k6,6rg -k2,2rg $data1 |\
                    awk '{print NR-1, $0}' > $data2
                plotOption1="
plot '$data2'  using 3:xtic(1) ls 1 t 'IDT',\
     ''       using 7:xtic(1) ls 6 t 'DIFF', \
    ''        using 6:xtic(1) ls 5 t 'INV_SHIFT', \
     ''       using 5:xtic(1) ls 4 t 'INV', \
     ''       using 4:xtic(1) ls 3 t 'SHIFT'
"
                ;;
            INV|inv)
                sort -k4,4rg -k3,3rg -k5,5rg -k6,6rg -k2,2rg $data1 |\
                    awk '{print NR-1, $0}' > $data2
                plotOption1="
plot '$data2'  using 3:xtic(1) ls 1 t 'IDT',\
     ''       using 7:xtic(1) ls 6 t 'DIFF', \
    ''        using 6:xtic(1) ls 5 t 'INV_SHIFT', \
     ''       using 4:xtic(1) ls 3 t 'SHIFT',\
     ''       using 5:xtic(1) ls 4 t 'INV'
"
                ;;
            INV_SHIFT|inv_shift)
                sort -k5,5rg -k4,4rg -k3,3rg -k6,6rg -k2,2rg $data1 |\
                    awk '{print NR-1, $0}' > $data2
                plotOption1="
plot '$data2'  using 3:xtic(1) ls 1 t 'IDT',\
     ''       using 7:xtic(1) ls 6 t 'DIFF', \
     ''       using 4:xtic(1) ls 3 t 'SHIFT',\
     ''       using 5:xtic(1) ls 4 t 'INV',\
    ''        using 6:xtic(1) ls 5 t 'INV_SHIFT'
"
                ;;
            DIFF|diff)
                sort -k6,6rg -k4,4rg -k3,3rg -k5,5rg -k2,2rg $data1 |\
                    awk '{print NR-1, $0}' > $data2
                plotOption1="
plot \
     '$data2'  using 3:xtic(1) ls 1 t 'IDT',\
     '$data2'  using 6:xtic(1) ls 5 t 'INV_SHIFT',\
     '$data2'  using 4:xtic(1) ls 3 t 'SHIFT',\
     '$data2'  using 5:xtic(1) ls 4 t 'INV',\
     '$data2'  using 7:xtic(1) ls 6 t 'DIFF'
"
                ;;
        esac
    fi
fi

outputSetting=
case $outputStyle in
    term*)
    outputSetting=
    ;;
    png)
    outfile=$outpath/${rootname}_${ordertype}_minnumseq${min_numseq}_minperidt${min_peridt}.png
    outputSetting="
set terminal png enhanced
set output '$outfile' 
    "
    ;;
    eps)
    outfile=$outpath/${rootname}_${ordertype}_minnumseq${min_numseq}_minperidt${min_peridt}.eps
    pdffile=$outpath/${rootname}_${ordertype}_minnumseq${min_numseq}_minperidt${min_peridt}.pdf
    outputSetting="
set term postscript eps color
set output '$outfile' 
    "
    ;;
esac

if [ -s $data2 ]; then 
/usr/bin/gnuplot -persist<<EOF 
$outputSetting
#set term postscript eps color
#set output "$outfile"

set style line 1 lt 1 lw 1 lc rgb "green"
set style line 3 lt 1 lw 1 lc rgb "violet"
set style line 4 lt 1 lw 1 lc rgb "blue" 
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "black"

set multiplot layout 2, 1 title ""

set lmargin at screen 0.12
set rmargin at screen 0.8
set bmargin at screen 0.4
set tmargin at screen 0.9

set ylabel "Composition in percentages"
set key outside top
#set key reverse
set style data histograms
set style histogram rowstacked
set boxwidth 1 relative
#set style fill solid 1.0 border -1
set style fill solid noborder
set xrange[0:]
set yrange [0:100]
unset xtics
#set size 1.0, 0.5
$plotOption1

#=============================
#    Plot 2
#=============================
set lmargin at screen 0.12
set rmargin at screen 0.8
set bmargin at screen 0.1 
set tmargin at screen 0.36

set xrange [0:$numPfam]
set autoscale y
set ytics nomirror
set xtics nomirror
set xtics
set xlabel "$numPfam Pfam families"
set ylabel "Number of sequences"
set style data linespoints
unset key
plot '$data2' using 1:2 ls 6 lt -1 pt 6 ps 0

#replot
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

# echo $plotOption1
#echo $tmpdir
rm -rf $tmpdir

