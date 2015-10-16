#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
#  #Idx  SeqIDT       IDT       INV    TM2GAP    TM2SEQ TM2GAP_AND_TM2SEQ   Maximum Occurrence
# 0       0-10    23.378     4.501    15.654    13.330    43.138    43.138      59702
# 1      10-20    50.891     4.041    19.435    12.237    13.396    50.891      23560
# 2      20-30    70.897     2.489    11.571    10.671     4.372    70.897      25753
# 3      30-40    76.495     2.398     7.739    10.697     2.670    76.495      18762
# 4      40-50    78.959     2.145     4.021    13.445     1.430    78.959      10071
# 5      50-60    83.021     2.340     2.340    11.119     1.181    83.021       4488
# 6      60-70    86.222     1.385     2.450     9.197     0.746    86.222       2816
# 7      70-80    89.414     0.806     1.666     7.523     0.591    89.414       1861
# 8      80-90    89.466     0.964     1.558     7.789     0.223    89.466       1348
# 9     90-100    95.541     0.398     0.318     3.583     0.159    95.541       1256

usage="
Usage: $0 datafile
Description: 
Options:
    -o, -outstyle STR   set the output, (default: eps)
                        can be png, term, eps
    -outpath      DIR   set the output path, default: the same as input file
    -numcase y|n        show/not show number of cases, default: not show
    -twoclass           show just identical/different
    -plot1              Print merged bin, 0-10,10-20,20-30,30-40,40-100
    -multiplot          Show number of cases as multiplot
    -h, --help          print this help message and exit

Created 2013-10-07, updated 2013-10-07, Nanjiang Shu 

Example:
    $0 datafile -twoclass -multiplot -plot1
    $0 datafile datafile -plot1
"
PrintHelp(){
    echo "$usage"
}

GetYTIC2() { #maxNum
    local maxnum=$1
    #local maxNumExpample=`awk '{print $NF}' $dataFile | sort -rg | head -n 1`
    local ytic2=
    ((ytic2=maxnum/4))
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
CreateDataFile1(){ #{{{
    local infile="$1"
    local outfile="$2"
awk -v method=$method -v idx=$idx '
BEGIN{
    totalsum = 0;
    COL=0;
    STARTINDEX = 4
    for (j=3;j<=9;j++){
        subsum[i] = 0;
    }
}

{
    if ($1>=STARTINDEX){
        COL=NF;
        totalsum += $NF;
        for (i=3;i<=NF-2;i++){
            cnt[i] = $i*$NF/100;
            subsum[i] +=  cnt[i]
#             print int(cnt[i]), int(subsum[i]);
        }
    }else if ($1 != "#sum"){
        print
    }
}
END{
    for (i=3;i<=COL-2;i++){
        freq[i] = subsum[i]/totalsum;
    }
    printf("%d ", STARTINDEX);
    printf("%d-100 ", STARTINDEX*10);
    maxf = -1.0
    for (i=3;i<=COL-2;i++){
        printf("%6.3f ", freq[i]*100);
        if (freq[i]*100 > maxf){
            maxf = freq[i]*100
        }
    }
    printf("%6.3f %d \n", maxf, totalsum);
}' $infile > $outfile
}
#}}}
CreateDataFile2(){ #{{{
    local infile="$1"
    local outfile="$2"
awk -v method=$method -v idx=$idx '
BEGIN{
#for merging  0-10 and 10-20
    totalsum0 = 0
    for (j=3;j<=9; j++){
        subsum0[i] = 0;
    }

    totalsum = 0;
    COL=-1;
    STARTINDEX = 4


    for (j=3;j<=9;j++){
        subsum[i] = 0;
    }
}

{
    if ($1 == 0 || $1 == 1){
        if(COL == -1){ COL=NF;}
        totalsum0 += $NF;
        for (i=3;i<=NF-2;i++){
            subsum0[i] +=  $i*$NF/100
        }
    }else if ($1 == 2){
        printf("%d ", 0);
        printf("0-20 ");
        maxf = -1.0
        for (i=3;i<=COL-2;i++){
            freq[i] = subsum0[i]/totalsum0;
            printf("%6.3f ", freq[i]*100);
            if (freq[i]*100 > maxf){
                maxf = freq[i]*100
            }
        }
        printf("%6.3f %d \n", maxf, totalsum0);

        printf("%d ", $1-1);
        printf("%s ", $2);
        for (i=3;i<=NF-2;i++){
            printf("%6.3f ", $i);
        }
        printf("%6.3f %d \n", $(NF-1), $NF);


    }else if ($1>=STARTINDEX){
        if(COL == -1){ COL=NF;}
        totalsum += $NF;
        for (i=3;i<=NF-2;i++){
            cnt[i] = $i*$NF/100;
            subsum[i] +=  cnt[i]
#             print int(cnt[i]), int(subsum[i]);
        }
    }else if (substr($1,1,1) != "#"){
        printf("%d ", $1-1);
        printf("%s ", $2);
        for (i=3;i<=NF-2;i++){
            printf("%6.3f ", $i);
        }
        printf("%6.3f %d \n", $(NF-1), $NF);
    }else if ($1 != "#sum"){
        print $0
    }
}
END{
    for (i=3;i<=COL-2;i++){
        freq[i] = subsum[i]/totalsum;
    }
    printf("%d ", STARTINDEX-1);
    printf("%d-100 ", STARTINDEX*10);
    maxf = -1.0
    for (i=3;i<=COL-2;i++){
        printf("%6.3f ", freq[i]*100);
        if (freq[i]*100 > maxf){
            maxf = freq[i]*100
        }
    }
    printf("%6.3f %d \n", maxf, totalsum);
}' $infile > $outfile
}
#}}}
MakePlot_shownumcase(){ #{{{  #showing number of cases
    local dataFile="$1"
    local outputSetting=
    local basename=`basename "$dataFile"`
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
    pdffile=$outpath/$basename.pdf
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

set title "$title"
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages of topology comparison classes"  font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
set label 1 at 3.5, 110
set label 1 "Number of cases"

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
set grid ytics
#set bmargin 5

plot \
    '$dataFile' using 3:xtic(2)   title 'Identical' ls 1 , \
    ''          using 4:xtic(2)   title 'Inverted' ls 2, \
    ''          using 5:xtic(2)   title 'TM2GAP' ls 3, \
    ''          using 6:xtic(2)   title 'TM2SEQ' ls 4, \
    ''          using 7:xtic(2)   title 'TM2GAP and TM2SEQ' ls 5, \
    ''          using (\$1):(\$3+\$4+\$5+\$6+\$7+3):9 title '' with labels
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
MakePlot_no(){ #{{{  #not show number of cases
    local dataFile="$1"
    local basename=`basename "$dataFile"` 
    local outputSetting=
    local b1=$basename
    local plotOption=""

    if [ $isTwoClass -eq 0 ];then
        plotOption="'$dataFile' using 3:xtic(2)   title 'Identical' ls 1 , \
        ''          using 4:xtic(2)   title 'Inverted' ls 2, \
        ''          using 5:xtic(2)   title 'TM2GAP' ls 3, \
        ''          using 6:xtic(2)   title 'TM2SEQ' ls 4, \
        ''          using 7:xtic(2)   title 'Both TM2GAP and TM2SEQ' ls 5
        "
    else
        plotOption="'$dataFile' using (\$4+\$5+\$6+\$7):xtic(2)   title 'Different' ls 11 , \
        ''          using 3:xtic(2)   title 'Identical' ls 12
        "
        b1=$b1.2cls
    fi

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
set style line 2 lt 1 lw 1 lc rgb "blue" 
set style line 3 lt 1 lw 1 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "violet"
set style line 5 lt 1 lw 1 lc rgb "cyan" 
set style line 6 lt 1 lw 1 lc rgb "orange"
set style line 7 lt 1 lw 1 lc rgb "grey40"
set style line 8 lt 1 lw 1 lc rgb "grey70"
set style line 11 lt 1 lw 1 lc rgb "black"
set style line 12 lt 1 lw 1 lc rgb "grey90"

set title "$title"
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages of topology comparison classes"  font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
set yrange[0:100]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
set grid ytics

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
MakePlot_no_multiplot(){ #{{{  #show also number of cases as multiplot
    local dataFile="$1"
    local basename=`basename "$dataFile"` 
    local outputSetting=
    local b1=$basename
    local plotOption1=""
    local plotOption2=""
    local NF=`awk '{print NF}' $dataFile | head -n 1` # number of fields
    local N_record=`grep "^[^#]" $dataFile | wc -l`
    b1=$b1.mtp


    local maxNumExpample=`awk '{print $NF}' $dataFile | sort -rg | head -n 1`
    ytic2=`GetYTIC2 $maxNumExpample`

    if [ $isTwoClass -eq 0 ];then
        plotOption1="\
        '$dataFile' using 3:xtic(2)   title 'Identical' ls 1, \
        ''          using 4:xtic(2)   title 'Inverted'  ls 2, \
        ''          using 5:xtic(2)   title 'TM2GAP'    ls 3, \
        ''          using 8:xtic(2)   title 'Mixed'     ls 6, \
        ''          using 6:xtic(2)   title 'TM2SEQ'    ls 4, \
        ''          using 7:xtic(2)   title 'TM2SP'     ls 5  \
        "
    else
        plotOption1="\
        '$dataFile' using (\$4+\$5+\$6+\$7+\$8):xtic(2) title 'Different topology' ls 11 , \
        ''          using 3:xtic(2)   title 'Identical topology' ls 12\
        "
        b1=$b1.2cls
    fi
    plotOption2="'$dataFile' using $NF:xtic(2) ti '' ls 21 "

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
set style line 1  lt 1 lw 1 lc rgb "red"
set style line 2  lt 1 lw 1 lc rgb "blue" 
set style line 3  lt 1 lw 1 lc rgb "green"
set style line 4  lt 1 lw 1 lc rgb "violet"
set style line 5  lt 1 lw 1 lc rgb "cyan" 
set style line 6  lt 1 lw 1 lc rgb "orange"
set style line 7  lt 1 lw 1 lc rgb "grey40"
set style line 8  lt 1 lw 1 lc rgb "grey70"
set style line 11 lt 1 lw 1 lc rgb "black"
set style line 12 lt 1 lw 1 lc rgb "grey90"
set style line 21 lt 2 lw 2 pt 4 ps 2 lc rgb "black"

set multiplot layout 2, 1 title ""

set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.35
set tmargin at screen 0.85


set title "$title"  font "Arial, 24"
unset xlabel
unset xtics
#set xlabel "$xlabel"  font "Arial, 20"
set ylabel "$ylabel"  font "Arial, 20"
#set key top center outside
#set key bottom center outside
#set key top right box
#set key bottom center outside
#set key at screen 0.9, 0.925 font ",18"
set key at screen 0.35, 0.96 font ",18"
set key spacing 1.3
set key invert
set autoscale y
set yrange[0:100]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics font ',18'
set grid ytics
set autoscale x
set xrange[-1:$N_record]

plot \
    $plotOption1

#=============================
#    Plot 2
#=============================
set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.10 
set tmargin at screen 0.30

set title ""
set xlabel "$xlabel"  font "Arial, 18"
set xtics font "Arial,20"
unset x2tics
set ylabel "$ylabel_plot2"
set ytics font 'Arial,18'
set grid noxtics ytics
set style data linespoints
#set style data histogram
#set style nofill solid border -1
#set style fill pattern border
set autoscale x
set autoscale y
set yrange[0:]
set ytics $ytic2
set xrange[-1:$N_record]
unset key

plot \
    $plotOption2

unset multiplot

EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            pdfcrop $pdffile
            echo "Histogram image output to $pdffile"
            ;;
        *) echo "Histogram image output to $outfile" ;;
    esac

}
#}}}
MakePlot_no_multiplot_noidt(){ #{{{  #show also number of cases as multiplot
    # do not show the grey bar of identical topology
    local dataFile="$1"
    local basename=`basename "$dataFile"` 
    local outputSetting=
    local b1=$basename.noidt
    local plotOption1=""
    local plotOption2=""
    local NF=`awk '{print NF}' $dataFile | head -n 1` # number of fields
    local N_record=`grep "^[^#]" $dataFile | wc -l`
    b1=$b1.mtp


    local maxNumExpample=`awk '{print $NF}' $dataFile | sort -rg | head -n 1`
    ytic2=`GetYTIC2 $maxNumExpample`

    if [ $isTwoClass -eq 0 ];then
        plotOption1="\
        '$dataFile' using 3:xtic(2)   title 'Identical' ls 1, \
        ''          using 4:xtic(2)   title 'Inverted'  ls 2, \
        ''          using 5:xtic(2)   title 'TM2GAP'    ls 3, \
        ''          using 8:xtic(2)   title 'Mixed'     ls 6, \
        ''          using 6:xtic(2)   title 'TM2SEQ'    ls 4, \
        ''          using 7:xtic(2)   title 'TM2SP'     ls 5  \
        "
    else
        plotOption1="\
        '$dataFile' using (\$4+\$5+\$6+\$7+\$8):xtic(2) title 'Different topology' ls 11\
        "
        b1=$b1.2cls
    fi
    plotOption2="'$dataFile' using $NF:xtic(2) ti '' ls 21 "

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
set style line 1  lt 1 lw 1 lc rgb "red"
set style line 2  lt 1 lw 1 lc rgb "blue" 
set style line 3  lt 1 lw 1 lc rgb "green"
set style line 4  lt 1 lw 1 lc rgb "violet"
set style line 5  lt 1 lw 1 lc rgb "cyan" 
set style line 6  lt 1 lw 1 lc rgb "orange"
set style line 7  lt 1 lw 1 lc rgb "grey40"
set style line 8  lt 1 lw 1 lc rgb "grey70"
set style line 11 lt 1 lw 1 lc rgb "black"
set style line 12 lt 1 lw 1 lc rgb "grey90"
set style line 21 lt 2 lw 2 pt 4 ps 2 lc rgb "black"

set multiplot layout 2, 1 title ""

set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.35
set tmargin at screen 0.85


set title "$title"  font "Arial, 24"
unset xlabel
unset xtics
#set xlabel "$xlabel"  font "Arial, 20"
set ylabel "$ylabel"  font "Arial, 20"
#set key top center outside
#set key bottom center outside
#set key top right box
#set key bottom center outside
#set key at screen 0.9, 0.925 font ",18"
set key at screen 0.35, 0.96 font ",18"
set key spacing 1.3
set key invert
set autoscale y
set yrange[0:65]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics font ',18'
set grid ytics
set autoscale x
set xrange[-1:$N_record]

plot \
    $plotOption1

#=============================
#    Plot 2
#=============================
set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.10 
set tmargin at screen 0.30

set title ""
set xlabel "$xlabel"  font "Arial, 18"
set xtics font "Arial,20"
unset x2tics
set ylabel "$ylabel_plot2"
set ytics font 'Arial,18'
set grid noxtics ytics
set style data linespoints
#set style data histogram
#set style nofill solid border -1
#set style fill pattern border
set autoscale x
set autoscale y
set yrange[0:]
set ytics $ytic2
set xrange[-1:$N_record]
unset key

plot \
    $plotOption2

unset multiplot

EOF

    case $outputStyle in
        eps) my_epstopdf $outfile
            pdfcrop $pdffile
            echo "Histogram image output to $pdffile"
            ;;
        *) echo "Histogram image output to $outfile" ;;
    esac

}
#}}}


MakePlot_no_invert_multiplot(){ #{{{  #show also number of cases as multiplot
    local dataFile="$1"
    local basename=`basename "$dataFile"` 
    local outputSetting=
    local b1=$basename
    local plotOption1=""
    local plotOption2=""
    local NF=`awk '{print NF}' $dataFile | head -n 1` # number of fields
#    local maxNumExpample=`awk '{print $NF}' $dataFile | sort -rg | head -n 1`
#    local ytic2=
#    ((ytic2=maxNumExpample/3))
#    if [ $ytic2 -eq 0 ];then
#        ytic2=1
#    fi
    local N_record=`grep "^[^#]" $dataFile | wc -l`
    b1=$b1.revert.mtp

    plotOption1="\
    '$dataFile' using 4:xtic(2)   title 'Inverted' $fs2 ls 2, \
    ''          using 5:xtic(2)   title 'TM2GAP'   $fs3 ls 3, \
    ''          using 8:xtic(2)   title 'Mixed'    $fs6 ls 6, \
    ''          using 6:xtic(2)   title 'TM2SEQ'   $fs4 ls 4, \
    ''          using 7:xtic(2)   title 'TM2SP'    $fs5 ls 5  \
    "
    plotOption2="'$dataFile' using $NF:xtic(2) ti '' ls 21 "

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
# set style line 1  lt 1 lw 1 lc rgb "red"
# set style line 2  lt 1 lw 1 lc rgb "blue" 
# set style line 3  lt 1 lw 1 lc rgb "green"
# set style line 4  lt 1 lw 1 lc rgb "violet"
# set style line 5  lt 1 lw 1 lc rgb "cyan" 
# set style line 6  lt 1 lw 1 lc rgb "orange"
# black and white
set style line 1  lt 1 lw 1 lc rgb "grey50"
set style line 2  lt 1 lw 1 lc rgb "grey50" 
set style line 3  lt 1 lw 1 lc rgb "grey50"
set style line 4  lt 1 lw 1 lc rgb "grey50"
set style line 5  lt 1 lw 1 lc rgb "grey50" 
set style line 6  lt 1 lw 1 lc rgb "grey50"

set style line 7  lt 1 lw 1 lc rgb "grey40"
set style line 8  lt 1 lw 1 lc rgb "grey70"
set style line 11 lt 1 lw 1 lc rgb "black"
set style line 12 lt 1 lw 1 lc rgb "grey90"


set style line 21 lt 2 lw 2 pt 4 ps 2 lc rgb "black"

set multiplot layout 2, 1 title ""

set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.30
set tmargin at screen 0.80


set title "$title"
unset xlabel
unset xtics
set ylabel "$ylabel"  font "Arial, 16"
#set key top center outside
#set key bottom center outside
set key top right boxed
set key invert
#set yrange[0:100]
set autoscale y

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set ytics font 'Times-Roman,16'
set grid ytics
set autoscale x
set xrange[-1:$N_record]

plot $plotOption1

#=============================
#    Plot 2
#=============================
set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.10 
set tmargin at screen 0.25

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set xtics
unset x2tics
set ylabel "$ylabel_plot2"
set ytics font 'Arial,16'
set grid noxtics ytics
set style data linespoints
#set style data histogram
#set style nofill solid border -1
#set style fill pattern border
set autoscale y
set autoscale y
set yrange[0:]
#set ytics $ytic2
set ytics 400
#set mytics 5
set autoscale x
set xrange[-1:$N_record]
unset key

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
MakePlot_no_invert(){ #{{{  #not show number of cases
    local dataFile="$1"
    local basename=`basename "$dataFile"` 
    basename=$basename.revert
    local outputSetting=""
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
    pdffile=$outpath/$basename.pdf
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

set title "$title"
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "Percentages of topology comparison classes"  font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
set yrange[0:32]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
set grid ytics

plot \
    '$dataFile' using 4:xtic(2)   title 'Inverted' ls 2, \
    ''          using 5:xtic(2)   title 'TM2GAP' ls 3, \
    ''          using 6:xtic(2)   title 'TM2SEQ' ls 4, \
    ''          using 7:xtic(2)   title 'TM2GAP and TM2SEQ' ls 5
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

title=
outputStyle=eps
outpath=
dataFile=
xlabel="Sequence identity"
ylabel="Relative frequency (%)"
ylabel_plot2="Number of pairs"
isNonOptionArg=false
isShowNumCase=0
isPlot1=0
isTwoClass=0
isMultiplot=0

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
            -plot1|--plot1) isPlot1=1;;
            -twoclass|--twoclass) isTwoClass=1;;
            -multiplot|--multiplot) isMultiplot=1;;
            -showcase|--showcase) topt=$2; 
                if [ "${topt:0:1}" == "y" -o  "${topt:0:1}" == "Y" ]; then
                    isShowNumCase=1
                else
                    isShowNumCase=0
                fi
                shift;
                ;;
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
dataFile2=$stemname.2.txt

if [ $isPlot1 -eq 1 ];then
    CreateDataFile1 "$dataFile" "$dataFile1"
    CreateDataFile2 "$dataFile" "$dataFile2"
fi


if [ $isShowNumCase -eq 1 ]; then
    MakePlot_shownumcase $dataFile
    if [ $isPlot1 -eq 1 ];then
        MakePlot_shownumcase $dataFile1
        MakePlot_shownumcase $dataFile2
    fi
else
    if [ $isMultiplot -eq 1 ];then
        MakePlot_no_multiplot $dataFile
        if [ $isTwoClass -eq 1 ];then
            # added 2014-10-09, for twoclass display, do not show identical
            # topology
            MakePlot_no_multiplot_noidt $dataFile
        fi
    else
        MakePlot_no $dataFile
    fi

    if [ $isPlot1 -eq 1 ];then
        for file in $dataFile1 $dataFile2; do
            MakePlot_no_invert $file
            if [ $isMultiplot -eq 1 ];then
                MakePlot_no_multiplot $file
                MakePlot_no_invert_multiplot $file
            else
                MakePlot_no $file
            fi
        done
    fi
fi

