#!/bin/bash
# Filename: plotGOFreqAna.sh 
# Description:
#   Draw histogram for GO frequency analysis

usage="
Usage: $0 datafiles -outname STR
Description:
Options:
    -outstyle png|term|eps  set the output, default = terminal
    -plotnormal
    -plotwithconf
    -plotwithpvalue
    -h|--help               print this help message and exit

Created 2013-09-12, updated 2013-09-18, Nanjiang Shu 
"

os=`uname -s`
case $os in 
    Darwin*)
        sed_exec=gsed;;
    *)
        sed_exec=sed;;
esac
echo "sed_exec=$sed_exec"

PrintHelp() {
    echo "$usage"
}

PlotWithGOID() { #{{{
    local outputSetting=
    local outfile=
    local pdffile=
    local basename=$basename.GOID
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$basename.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$basename.eps
        pdffile=$basename.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac

    #trap 'rm -rf "$tmpdir"' INT TERM EXIT


    paste $filelist | awk -F"\t" -v NUM=$numFile '{for(j=0;j<NUM;j++){printf "%s\t%f\t", $(4*j+1), $(4*j+4) }  printf "\n" }' | awk -F "\t" -v NUM=$numFile '{isprint=0; for(j=0;j<NUM;j++){if($(j*2+2)>0.01){isprint=1}}; if(isprint==1){print}}' > $basename.txt
#cat $basename.txt

    local plotOption=
    local title=
    local i1=
    local i2=
    local ils=
    for ((i=0;i<numFile;i++));do
        title=${titleList[$i]}
        ((i1=$i*2+1))
        ((i2=$i*2+2))
        case $title in 
            *IDT*) ils=1; title="Identical";;
            *INV*) ils=2; title="Inverted";;
            *TM2GAP*) ils=3; title="TM2GAP";;
            *TM2SEQ*) ils=4; title="TM2SEQ";;
            *SP2TM*) ils=5; title="SP2TM";;
            *TM2SP*) ils=5; title="TM2SP";;
            *Other*|other*) ils=6; title=Other;;
            *all*) ils=7; title=All;;
        esac
        if [ $i -eq 0 ];then
            plotOption="'$basename.txt' using $i2:xtic($i1) ti '$title' ls $ils"
        else
            plotOption="$plotOption, '' using $i2:xtic($i1) ti '$title' ls $ils"
        fi
    done
    echo "$plotOption"

/usr/bin/gnuplot -persist<<EOF 
$outputSetting

set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "pink"
set style line 6 lt 1 lw 1 lc rgb "cyan"   #all
set style line 7 lt 1 lw 1 lc rgb "grey70"

set style line 8 lt 1 lw 1 lc rgb "grey40"
set style line 9 lt 1 lw 1 lc rgb "grey70"
#set datafile separator "\t"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"

set style data histogram 
set style histogram clustered gap 1
set style fill solid border -1
set boxwidth 1
set tics font 'Times-Roman,16'
set xtics rotate by -30

plot \
    $plotOption
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
PlotWithGOTerm() { #{{{
    local outputSetting=
    local outfile=
    local pdffile=
    local basename=$basename.GOTERM
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$basename.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$basename.eps
        pdffile=$basename.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac

    #trap 'rm -rf "$tmpdir"' INT TERM EXIT


    #paste $filelist | $sed_exec 's/ *\t */\t/g' |sed 's/activity//g' | awk -F"\t" -v NUM=$numFile '{for(j=0;j<NUM;j++){printf "%s\t%f\t", $(4*j+2), $(4*j+4) }  printf "\n" }' | awk -F "\t" -v NUM=$numFile '{isprint=0; for(j=0;j<NUM;j++){if($(j*2+2)>0.01){isprint=1}}; if(isprint==1){print}}' > $basename.txt
    paste $filelist | $sed_exec 's/ *\t */\t/g' | awk -F"\t" -v NUM=$numFile '{for(j=0;j<NUM;j++){printf "%s\t%f\t", $(4*j+2), $(4*j+4) }  printf "\n" }' | awk -F "\t" -v NUM=$numFile '{isprint=0; for(j=0;j<NUM;j++){if($(j*2+2)>0.01){isprint=1}}; if(isprint==1){print}}' > $basename.txt
#cat $basename.txt

    local plotOption=
    local title=
    local i1=
    local i2=
    local ils=
    for ((i=0;i<numFile;i++));do
        title=${titleList[$i]}
        ((i1=$i*2+1))
        ((i2=$i*2+2))
        case $title in 
            IDT*) ils=1; title="Identical";;
            INV*) ils=2; title="INV|noSP";;
            TM2GAP*) ils=3; title="TM2GAP|noSP";;
            TM2SEQ*) ils=4; title="TM2SEQ|noSP";;
            SP2TM*) ils=5; title="SP2TM";;
            SP2TM*) ils=5; title=SP2TM;;
            Other*|other*) ils=6; title=Other;;
            all*) ils=7; title=All;;
        esac
        if [ $i -eq 0 ];then
            plotOption="'$basename.txt' using $i2:xtic($i1) ti '$title' ls $ils"
        else
            plotOption="$plotOption, '' using $i2:xtic($i1) ti '$title' ls $ils"
        fi
    done
    echo "$plotOption"

/usr/bin/gnuplot -persist<<EOF 
$outputSetting

set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "pink"
set style line 6 lt 1 lw 1 lc rgb "cyan"   #all
set style line 7 lt 1 lw 1 lc rgb "grey70"

set style line 8 lt 1 lw 1 lc rgb "grey40"
set style line 9 lt 1 lw 1 lc rgb "grey70"
set datafile separator "\t"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"

set style data histogram 
set style histogram clustered gap 1
set style fill solid border -1
set boxwidth 1
set tics font 'Times-Roman,16'
set xtics rotate by -30 

plot \
    $plotOption
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
PlotByTopo() { #{{{
    local file=$1
    local basename=`basename $file`
    local outputSetting=
    local outfile=
    local pdffile=
    local title=${file%.*}
    title=${title%.*}
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$basename.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$basename.eps
        pdffile=$basename.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac

    local tmpfile=$(mktemp /tmp/tmp.plotGOFreqAna.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    trap 'rm -f "$tmpfile"' INT TERM EXIT

    #cat $file | $sed_exec 's/ *\t */\t/g' |sed 's/activity//g' |awk -F "\t"  '{if($4>=0.01){print}}'   > $tmpfile
    cat $file | $sed_exec 's/ *\t */\t/g' | awk -F "\t"  '{if($4>=0.01){print}}'   > $tmpfile
    cat $tmpfile

    local ils=
    case $title in 
        IDT*) ils=1; title="Identical";;
        INV*) ils=2; title="INV|noSP";;
        TM2GAP*) ils=3; title="TM2GAP|noSP";;
        TM2SEQ*) ils=4; title="TM2SEQ|noSP";;
        SP2TM*) ils=5; title="SP2TM";;
        SP2TM*) ils=5; title=SP2TM;;
        Other*|other*) ils=6; title=Other;;
        all*) ils=7; title=All;;
    esac
    local plotOption="'$tmpfile' using 4:xtic(2) title '$title' ls $ils"


/usr/bin/gnuplot -persist<<EOF 
$outputSetting

set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "pink"
set style line 6 lt 1 lw 1 lc rgb "cyan"   #all
set style line 7 lt 1 lw 1 lc rgb "grey70"

set style line 8 lt 1 lw 1 lc rgb "grey40"
set style line 9 lt 1 lw 1 lc rgb "grey70"
set datafile separator "\t"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"

set style data histogram 
set style histogram clustered gap 1
set style fill solid border -1
set boxwidth 1
set tics font 'Times-Roman,16'
set xtics rotate by -30 

plot \
    $plotOption
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
}
#}}}
PlotByTopoWithAll() { #{{{
    local file=$1
    local basename=`basename $file`
    basename=$basename.withall
    local outputSetting=
    local outfile=
    local pdffile=
    local title=${file%.*}
    title=${title%.*}
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$basename.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$basename.eps
        pdffile=$basename.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac

    local tmpfile=$(mktemp /tmp/tmp.plotGOFreqAna.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    trap 'rm -f "$tmpfile"' INT TERM EXIT

    #cat $file | $sed_exec 's/ *\t */\t/g' |sed 's/activity//g' |awk -F "\t"  '{if($4>=0.01){print}}'   > $tmpfile
    cat $file | $sed_exec 's/ *\t */\t/g' | awk -F "\t"  '{if($4>=0.01){print}}'   > $tmpfile
    cat $tmpfile

    local ils=
    case $title in 
        IDT*) ils=1; title="Identical";;
        INV*) ils=2; title="INV|noSP";;
        TM2GAP*) ils=3; title="TM2GAP|noSP";;
        TM2SEQ*) ils=4; title="TM2SEQ|noSP";;
        SP2TM*) ils=5; title="SP2TM";;
        SP2TM*) ils=5; title=SP2TM;;
        Other*|other*) ils=6; title=Other;;
        all*) ils=7; title=All;;
    esac
    local plotOption="'$tmpfile' using 4:xtic(2) title '$title' ls $ils, '$allfile' using 4:xtic(2) title 'All' ls 7"
    echo "$plotOption"


/usr/bin/gnuplot -persist<<EOF 
$outputSetting

set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "pink"
set style line 6 lt 1 lw 1 lc rgb "cyan"   #all
set style line 7 lt 1 lw 1 lc rgb "grey70"

set style line 8 lt 1 lw 1 lc rgb "grey40"
set style line 9 lt 1 lw 1 lc rgb "grey70"
set datafile separator "\t"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"

set style data histogram 
set style histogram clustered gap 1
set style fill solid border -1
set boxwidth 1
set tics font 'Times-Roman,16'
set xtics rotate by -30 

plot \
    $plotOption
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
}
#}}}
PlotByTopoMinusAll() { #{{{
    local file=$1
    local basename=`basename $file`
    basename=$basename.minusall
    local outputSetting=
    local outfile=
    local pdffile=
    local title=${file%.*}
    title=${title%.*}
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$basename.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$basename.eps
        pdffile=$basename.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac

    local tmpfile=$(mktemp /tmp/tmp.plotGOFreqAna.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    local tmpfile1=$(mktemp /tmp/tmp.plotGOFreqAna.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }   
    trap 'rm -f "$tmpfile"' INT TERM EXIT
    trap 'rm -f "$tmpfile1"' INT TERM EXIT
    awk -F "\t"  '{if($4>=0.01){print}}'  $file > $tmpfile1
    paste $tmpfile1 $allfile  |awk -F "\t"  '{if($4>=0.01){print}}'  | awk -F"\t" '{print $1 "\t" $2 "\t"  $3 "\t"  $4-$8 "\t" $4 "\t" $8}' > $tmpfile
    cat $tmpfile

    local ils=
    case $title in 
        IDT*) ils=1; title="Identical";;
        INV*) ils=2; title="INV|noSP";;
        TM2GAP*) ils=3; title="TM2GAP|noSP";;
        TM2SEQ*) ils=4; title="TM2SEQ|noSP";;
        SP2TM*) ils=5; title="SP2TM";;
        SP2TM*) ils=5; title=SP2TM;;
        Other*|other*) ils=6; title=Other;;
        all*) ils=7; title=All;;
    esac
    local plotOption="'$tmpfile' using 4:xtic(2) title '$title' ls $ils"
    echo "$plotOption"


/usr/bin/gnuplot -persist<<EOF 
$outputSetting

set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "pink"
set style line 6 lt 1 lw 1 lc rgb "cyan"   #all
set style line 7 lt 1 lw 1 lc rgb "grey70"

set style line 8 lt 1 lw 1 lc rgb "grey40"
set style line 9 lt 1 lw 1 lc rgb "grey70"
set datafile separator "\t"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"

set style data histogram 
set style histogram clustered gap 1
set style fill solid border -1
set boxwidth 1
set tics font 'Times-Roman,16'
set xtics rotate by -30 

plot \
    $plotOption
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
}
#}}}
PlotNormal() { #{{{
    local file=$1
    local basename=`basename $file`
    basename=$basename.basic
    local outputSetting=
    local outfile=
    local pdffile=
    local title=${file%.*}
    title=${title%.*}
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$basename.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$basename.eps
        pdffile=$basename.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac

    local tmpfile=$(mktemp /tmp/tmp.plotGOFreqAna.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    trap 'rm -f "$tmpfile"' INT TERM EXIT
    #cat $file | $sed_exec 's/ *\t */\t/g' |sed 's/activity//g' > $tmpfile
    cat $file | $sed_exec 's/ *\t */\t/g' > $tmpfile
    cat $tmpfile

    local ils=
    case $title in 
        */*IDT*) ils=1; title="Identical";;
        */*INV*) ils=2; title="INV|noSP";;
        */*TM2GAP*) ils=3; title="TM2GAP|noSP";;
        */*TM2SEQ*) ils=4; title="TM2SEQ|noSP";;
        */*SP2TM*) ils=5; title="SP2TM";;
        */*Mixed*|*/other*) ils=6; title=Mixed;;
        */*All*|*all*) ils=7; title=All;;
    *) ils=7;;
    esac
    local plotOption="'$tmpfile' using 4:xtic(2) w boxes title '$title' ls $ils"
    echo "$plotOption"


/usr/bin/gnuplot -persist<<EOF 
$outputSetting

set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "pink"
set style line 6 lt 1 lw 1 lc rgb "cyan"   #all
set style line 7 lt 1 lw 1 lc rgb "grey70"

set style line 8 lt 1 lw 1 lc rgb "grey40"
set style line 9 lt 1 lw 1 lc rgb "grey70"
set datafile separator "\t"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"

#set style data histogram 
#set style histogram clustered gap 1
set yrange[0:]
set style fill solid border -1
set boxwidth 1
set tics font 'Times-Roman,16'
set xtics rotate by -30 

plot \
    $plotOption
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
}
#}}}
PlotWithConf() { #{{{
    local file=$1
    local basename=`basename $file`
    basename=$basename.withconf
    local outputSetting=
    local outfile=
    local pdffile=
    local title=${file%.*}
    title=${title%.*}
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$basename.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$basename.eps
        pdffile=$basename.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac

    local tmpfile=$(mktemp /tmp/tmp.plotGOFreqAna.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    trap 'rm -f "$tmpfile"' INT TERM EXIT
    #cat $file | $sed_exec 's/ *\t */\t/g' |sed 's/activity//g' > $tmpfile
    cat $file | $sed_exec 's/ *\t */\t/g' > $tmpfile
#    cat $tmpfile

    local ils=
    case $title in 
        *DiffAll*) ils=5; title="Different";;
        *DiffWithout*) ils=5; title="Different (excluding TM2SP and TM2SEQ)";;
        *IDT*) ils=1; title="Identical";;
        *INV*) ils=2; title="Inverted";;
        *TM2GAP*) ils=3; title="TM2GAP";;
        *TM2SEQ*) ils=4; title="TM2SEQ";;
        *SP2TM*) ils=5; title="SP2TM";;
        *TM2SP*) ils=5; title="TM2SP";;
        *Mixed*|*/other*) ils=6; title=Mixed;;
        *All*|*all*) ils=7; title=All;;
        *) ils=7;;
    esac
    ils=7
    local plotOption="'$tmpfile' using 4:5:6:xtic(2)  title '$title' ls $ils"
    echo "$plotOption"


/usr/bin/gnuplot -persist<<EOF 
$outputSetting

set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "pink"
set style line 6 lt 1 lw 1 lc rgb "cyan"   #all
set style line 7 lt 1 lw 1 lc rgb "grey70"

set style line 8 lt 1 lw 1 lc rgb "grey40"
set style line 9 lt 1 lw 1 lc rgb "grey70"
set datafile separator "\t"

set grid y
#set key outside

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"

set style data histogram 
#set style histogram clustered gap -1
set style histogram errorbars gap 1 lw 1
set yrange[0:]
#set xrange[0:]
set style fill solid border -1
set boxwidth 1 relative
set tics font 'Times-Roman,16'
set xtics rotate by -30 

plot \
    $plotOption
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
}
#}}}
PlotWithPvalue() { #{{{
    local file=$1
    local basename=`basename $file`
    basename=$basename.wp
    local outputSetting=
    local outfile=
    local pdffile=
    local title=${file%.*}
    title=${title%.*}
    case $outputStyle in
        term*)
        outputSetting=
        ;;
        png)
        outfile=$basename.png
        outputSetting="
    set terminal png enhanced
    set output '$outfile' 
        "
        ;;
        eps)
        outfile=$basename.eps
        pdffile=$basename.pdf
        outputSetting="
    set term postscript eps color
    set output '$outfile' 
        "
        ;;
    esac

    local tmpfile=$(mktemp /tmp/tmp.plotGOFreqAna.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    trap 'rm -f "$tmpfile"' INT TERM EXIT
    #cat $file | awk -F"\t" '{printf("%d\t%s\t%s\\n(%s)\t%6.6g\t%6.6g\t%4.1e\n", NR-1, $1, $2, $1, $5, $8, $9) }' | $sed_exec 's/ *\t */\t/g' | $sed_exec 's/activity//g' > $tmpfile
    cat $file | awk -F"\t" '{printf("%d\t%s\t%s\\n(%s)\t%6.6g\t%6.6g\t%4.1e\n", NR-1, $1, $2, $1, $5, $8, $9) }' | $sed_exec 's/ *\t */\t/g'  > $tmpfile
    local avg=`awk -F "\t" '{if(NR==1) print $5}' $tmpfile`
    local yshift=`awk -F"\t" 'BEGIN{max=-1}{if($4>max)max=$4}END{printf("%g", max/25.0)}' $tmpfile`
    #cat $tmpfile

    local ils=
    case $title in 
        *DiffAll*) ils=5; title="Different";;
        *DiffWithout*) ils=5; title="Different (excluding TM2SP and TM2SEQ)";;
        *IDT*) ils=1; title="Identical";;
        *INV*) ils=2; title="Inverted";;
        *TM2GAP*) ils=3; title="TM2GAP";;
        *TM2SEQ*) ils=4; title="TM2SEQ";;
        *SP2TM*) ils=5; title="SP2TM";;
        *TM2SP*) ils=5; title="TM2SP";;
        *DUP*) ils=5; title="Duplicated";;
        *Mixed*|*/other*) ils=6; title=Mixed;;
        *All*|*all*) ils=7; title=All;;
        *) ils=7;;
    esac
    ils=7
    local plotOption="'$tmpfile' using 4:xtic(3)  title '$title' ls $ils"
    local opt=
    opt="'' using (\$1):(\$4+$yshift):(\$6 < 0.05 ? gprintf(\"P < %3.0e\", \$6) : '') title '' with labels"
    plotOption="$plotOption, $opt"
    opt="$avg ls 11 title '' "
    plotOption="$plotOption, $opt"
    #echo "$plotOption"


/usr/bin/gnuplot -persist<<EOF 
$outputSetting

set style line 1 lt 1 lw 1 lc rgb "red"
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "pink"
set style line 6 lt 1 lw 1 lc rgb "cyan"   #all
set style line 7 lt 1 lw 1 lc rgb "grey70"

set style line 8 lt 1 lw 1 lc rgb "grey40"
set style line 9 lt 1 lw 1 lc rgb "grey70"

set style line 11 lt 2 lw 4 lc rgb "black"
set datafile separator "\t"
#set key outside

#set grid y

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel"  font "Arial, 16"

set style data histogram 
set style histogram clustered gap 1
#set style histogram errorbars gap 0 lw 1
set yrange[0:]
#set xrange[0:]
set style fill solid border -1
set boxwidth 1 relative
set tics font 'Times-Roman,16'
set xtics rotate by -30 

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
    rm -f $tmpfile
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outputStyle=eps
outname=
dataFileList=()
isNonOptionArg=false
xlabel="Molecular function"
ylabel="Normalized frequency"
isPlotNormal=0
isPlotWithConf=0
isPlotWithPvalue=0

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
        -plotnormal|--plotnormal)isPlotNormal=1;;
        -plotwithconf|--plotwithconf)isPlotWithConf=1;;
        -plotwithpvalue|--plotwithpvalue)isPlotWithPvalue=1;;
            -outname|--outname) outname=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFileList+=($1)
    fi
    shift
done

numFile=${#dataFileList[@]}

if [ $numFile -eq 0  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi

basename=$outname

titleList=()
filelist=""
file=
title=
for ((i=0;i<numFile;i++));do
    file=${dataFileList[$i]}
    title=${file%.*}
    title=${title%.*}
    titleList+=($title)
    filelist="$filelist $file"
done

if [ $isPlotWithConf -eq 1 ];then 
    for ((i=0;i<numFile;i++));do
        file=${dataFileList[$i]}
        PlotWithConf "$file"
    done
elif [ $isPlotWithPvalue -eq 1 ];then
    for ((i=0;i<numFile;i++));do
        file=${dataFileList[$i]}
        PlotWithPvalue "$file"
    done
elif [ $isPlotNormal -eq 1 ];then
    for ((i=0;i<numFile;i++));do
        file=${dataFileList[$i]}
        PlotNormal "$file"
    done
else
    if [ "$outname"  == "" ]; then 
        echo "Error! outname not set!"
        exit 1
    fi
    PlotWithGOID
    PlotWithGOTerm

    allfile=$(mktemp /tmp/tmp.plotGOFreqAna.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  
    trap 'rm -f "$allfile"' INT TERM EXIT
    for ((i=0;i<numFile;i++));do
        file=${dataFileList[$i]}
        case $file in 
            all*)
            #cat $file | $sed_exec 's/ *\t */\t/g' |sed 's/activity//g' |awk -F "\t"  '{if($4>=0.01){print}}'   > $allfile
            cat $file | $sed_exec 's/ *\t */\t/g' |awk -F "\t"  '{if($4>=0.01){print}}'   > $allfile
        esac
    done

    for ((i=0;i<numFile;i++));do
        file=${dataFileList[$i]}
        PlotByTopo "$file"
    done
    for ((i=0;i<numFile;i++));do
        file=${dataFileList[$i]}
        PlotByTopoWithAll "$file"
    done

    for ((i=0;i<numFile;i++));do
        file=${dataFileList[$i]}
        PlotByTopoMinusAll "$file"
    done

    rm -f $allfile
fi
