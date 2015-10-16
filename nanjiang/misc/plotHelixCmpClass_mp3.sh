#!/bin/bash
# draw the histogram given the file in the format
# #topology comparison of four states
# #
# format of the file to make the plot
#  #Idx  SeqIDT         0         1         2   Maximum Occurrence
# 0       0-10    69.431    17.082    13.488    69.431     503959
# 1      10-20    82.700    10.024     7.276    82.700     484708
# 2      20-30    90.896     5.825     3.279    90.896     262787
# 3      30-40    94.371     3.555     2.074    94.371     303315
# 4      40-50    95.357     2.809     1.833    95.357     240253
# 5      50-60    95.549     2.710     1.741    95.549      86855
# 6      60-70    96.863     1.500     1.637    96.863      48070
# 7      70-80    97.592     1.157     1.250    97.592      29029
# 8      80-90    97.947     0.969     1.083    97.947      18464
# 9     90-100    98.337     0.725     0.937    98.337      16960
# #sum            85.119     8.585     6.296    85.119    1994400   
usage="
Usage: $0 datafile
Description: 
Options:
    -outstyle,-o STR    Set the output, (default: eps)
                        can be png, term, eps
    -outpath DIR        Set the output path, (default: the same as datafile)
    -xlabel STR         Set xlabel, (default: Sequence identity)
    -ylabel STR         Set Ylabel
    -plot1              Print merged bin, 0-10,10-20,20-30,30-40,40-100
    -color              Show plots in color, (default: b/w)
    -multiplot          Print also the number of cases
    -h,--help           Print this help message and exit

Created 2013-10-07, updated 2014-10-09, Nanjiang Shu 
"
cmd="$*"
PrintHelp(){
    echo "$usage"
}
GetYTIC2() { #maxNum#{{{
    local maxnum=$1
    #local maxNumExpample=`awk '{print $NF}' $dataFile | sort -rg | head -n 1`
    local ytic2=
    ((ytic2=maxnum/3))
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
#}}}
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
CreateDataFile3(){ #{{{
    local infile="$1"
    local outfile="$2"
awk -v method=$method -v idx=$idx '
BEGIN{
    COL=-1;
    STARTINDEX = 4

#for merging  0-10 and 10-20
    totalsum0 = 0
    for (j=3;j<=9; j++){
        subsum0[i] = 0;
    }

#for merging  40-100
    totalsum = 0;
    for (j=3;j<=9;j++){
        subsum[i] = 0;
    }
#for merging  0-100
    totalsumAll = 0;
    for (j=3;j<=9;j++){
        subsumAll[i] = 0;
    }


}

{
    if (substr($1,1,1) != "#"){
        if(COL == -1){ COL=NF;}
        totalsumAll += $NF;
        for (i=3;i<=NF-2;i++){
            subsumAll[i] +=  $i*$NF/100
        }
    }
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

    # output for all proteins
    for (i=3;i<=COL-2;i++){
        freq[i] = subsumAll[i]/totalsumAll;
    }
    printf("%d ", STARTINDEX);
    printf("%d-100(All) ", 0);
    maxf = -1.0
    for (i=3;i<=COL-2;i++){
        printf("%6.3f ", freq[i]*100);
        if (freq[i]*100 > maxf){
            maxf = freq[i]*100
        }
    }
    printf("%6.3f %d \n", maxf, totalsumAll);
}' $infile > $outfile
}
#}}}
CreateDataFile4(){ #{{{
    ##0-20, 20-30, 30-40, 40-60, 60-100
    local infile="$1"
    local outfile="$2"
awk -v method=$method -v idx=$idx '
BEGIN{
#for merging  0-10 and 10-20
    totalsum0 = 0
    for (j=3;j<=9; j++){
        subsum0[i] = 0;
    }

#for merging  40-50 and 50-60
    totalsum1 = 0
    for (j=3;j<=9; j++){
        subsum1[i] = 0;
    }

    totalsum = 0;
    COL=-1;
    STARTINDEX = 6
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
    }else if ($1 == 4 || $1 == 5){
        if(COL == -1){ COL=NF;}
        totalsum1 += $NF;
        for (i=3;i<=NF-2;i++){
            subsum1[i] +=  $i*$NF/100
        }
    }else if ($1>=STARTINDEX){
        if(COL == -1){ COL=NF;}
        totalsum += $NF;
        for (i=3;i<=NF-2;i++){
            cnt[i] = $i*$NF/100;
            subsum[i] +=  cnt[i]
#             print int(cnt[i]), int(subsum[i]);
        }
    }


    if ($1 == 2){
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
    }else if ($1 == 6){
        printf("%d ", 0);
        printf("40-60 ");
        maxf = -1.0
        for (i=3;i<=COL-2;i++){
            freq[i] = subsum1[i]/totalsum1;
            printf("%6.3f ", freq[i]*100);
            if (freq[i]*100 > maxf){
                maxf = freq[i]*100
            }
        }
        printf("%6.3f %d \n", maxf, totalsum1);
    }else if ($1 == 3){
        printf("%d ", $1-1);
        printf("%s ", $2);
        for (i=3;i<=NF-2;i++){
            printf("%6.3f ", $i);
        }
        printf("%6.3f %d \n", $(NF-1), $NF);
    }else if (substr($1,1,1) == "#"){
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
MakePlot_shownumcase() {  #{{{
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

    /usr/bin/gnuplot -persist<<EOF 
$outputSetting
set style line 1 lt 1 lw 1 lc rgb "orange"
set style line 2 lt 1 lw 1 lc rgb "grey20"
set style line 3 lt 1 lw 1 lc rgb "purple"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
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
    '$dataFile' using 3:xtic(2)   title 'TM2TM' ls 1 , \
    ''          using 4:xtic(2)   title 'TM2GAP' ls 2, \
    ''          using 5:xtic(2)   title 'TM2SEQ' ls 3, \
    ''          using (\$1):(\$3+\$4+\$5+3):7 title '' with labels
EOF
    local pdffile=$outpath/$basename.pdf
    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Figure output to $pdffile"
#             if [ -s "$pdffile" ]; then
#                 rm -f "$outfile"
#             fi
            ;;
        *) echo "Figure output to $outfile" ;;
    esac
}
#}}}
MakePlot_no() {  #{{{
    local dataFile="$1"
    local basename=`basename "$dataFile"`
    local outfile=
    local pdffile=
    local outputSetting=""
    local plotOption=
    plotOption="\
    '$dataFile' using 3:xtic(2)   title 'TM2TM'  $fs1 ls 1, \
    ''          using 4:xtic(2)   title 'TM2GAP' $fs2 ls 2, \
    ''          using 5:xtic(2)   title 'TM2SEQ' $fs3 ls 3, \
    ''          using 6:xtic(2)   title 'TM2SP'  $fs4 ls 4  \
    "

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
set style line 1 lt 1 lw 1 lc rgb "grey30"
set style line 2 lt 1 lw 1 lc rgb "grey70"
set style line 3 lt 1 lw 1 lc rgb "grey70"
set style line 4 lt 1 lw 1 lc rgb "grey80"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
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

plot $plotOption
EOF
    local pdffile=$outpath/$basename.pdf

    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Figure output to $pdffile"
#             if [ -s "$pdffile" ]; then
#                 rm -f "$outfile"
#             fi
            ;;
        *) echo "Figure output to $outfile" ;;
    esac

}
#}}}
MakePlot_no_invert() {  #{{{
    local dataFile="$1"
    local basename=`basename "$dataFile"`
    local outfile=
    local pdffile=
    local outputSetting=""
    basename=$basename.revert
    local plotOption=
    plotOption="'$dataFile' using 4:xtic(2)   title 'TM2GAP' ls 2, \
    ''          using 5:xtic(2)   title 'TM2SEQ' ls 3, \
    ''          using 6:xtic(2)   title 'TM2SP' ls 4
    "
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
set style line 1 lt 1 lw 1 lc rgb "orange"
set style line 2 lt 1 lw 1 lc rgb "grey10"
set style line 3 lt 1 lw 1 lc rgb "grey30"
set style line 4 lt 1 lw 1 lc rgb "grey90"

set title ""
set xlabel "$xlabel"  font "Arial, 16"
set ylabel "$ylabel" font "Arial, 16"
set key top center outside
#set key bottom center outside
set key invert
set yrange[0:10]

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
set grid ytics

plot $plotOption
EOF
    local pdffile=$outpath/$basename.pdf

    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Figure output to $pdffile"
#             if [ -s "$pdffile" ]; then
#                 rm -f "$outfile"
#             fi
            ;;
        *) echo "Figure output to $outfile" ;;
    esac

}
#}}}
MakePlot_no_invert_multiplot() {  #{{{
    local dataFile="$1"
    local basename=`basename "$dataFile"`
    local outfile=
    local pdffile=
    local outputSetting=""
    local NF=`awk '{print NF}' $dataFile | head -n 1` # number of fields
    local N_record=`grep "^[^#]" $dataFile | wc -l`
    local plotOption1=
    local plotOption2=


    local maxNumExpample=`awk '/^[^#]/{print $NF}' $dataFile | sort -rg | head -n 1`
    ytic2=`GetYTIC2 $maxNumExpample`

    basename=${basename}.revert.mtp

    if [ $isColor -eq 1 ];then
        basename=${basename}.color
    fi

    if [ $isColor -eq 0 ];then
        plotOption1="\
        '$dataFile' using 4:xtic(2)   title 'TM2GAP' $fs2 ls 2, \
        ''          using 5:xtic(2)   title 'TM2SEQ' $fs3 ls 3, \
        ''          using 6:xtic(2)   title 'TM2SP'  $fs4 ls 4  \
        "
    else # in color
        plotOption1="\
        '$dataFile' using 4:xtic(2)   title 'TM2GAP' $fs1 ls 32, \
        ''          using 5:xtic(2)   title 'TM2SEQ' $fs1 ls 33, \
        ''          using 6:xtic(2)   title 'TM2SP'  $fs1 ls 34  \
        "
    fi
    plotOption2="'$dataFile' using $NF:xtic(2) ti '' ls 21"


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
set style line 1 lt 1 lw 1 lc rgb "grey50"
set style line 2 lt 1 lw 1 lc rgb "grey50"
set style line 3 lt 1 lw 1 lc rgb "grey50"
set style line 4 lt 1 lw 1 lc rgb "grey50"

set style line 32  lt 1 lw 1 lc rgb "navy"   # Helix TM2GAP
set style line 33  lt 1 lw 1 lc rgb "skyblue"    # Helix TM2SEQ
set style line 34  lt 1 lw 1 lc rgb "grey80"  # Helix TM2SP

#set style line 21 lt 1 lw 1 pt 4 lc rgb "black"
set style line 21 lt 2 lw 2 pt 4 ps 2 lc rgb "black"

set multiplot layout 2, 1 title ""

set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.28
set tmargin at screen 0.78

set title "$title" font "Arial, 20"
unset xlabel
unset xtics
set ylabel "$ylabel"  font "Arial, 16"
#set key top center outside
set key at screen 0.9, 0.96
set key invert
#set yrange[0:10]
set autoscale y

set style data histogram 
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set tics font 'Times-Roman,16'
set grid ytics
set ytics 1

unset xtics
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
set ylabel "Number of pairs"
set ytics font 'Arial,14'
set grid noxtics ytics
set style data linespoints
set autoscale y
set yrange[0:]
set ytics $ytic2
set autoscale x
set xrange[-1:$N_record]
unset key

plot $plotOption2
EOF

    local pdffile=$outpath/$basename.pdf
    case $outputStyle in
        eps) my_epstopdf $outfile
            echo "Figure output to $pdffile"
            if [ -s "$pdffile" ]; then
                pdfcrop $pdffile
#                 rm -f "$outfile"
            fi
            ;;
        *) echo "Figure output to $outfile" ;;
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
xlabel="Sequence identity"
ylabel="% of TM helix mapping groups"
title=
isColor=0

isNonOptionArg=false
isShowNumCase=0
isPlot1=0
isMultiplot=0

fs1="fs solid"
fs2="fs solid 0.33 "
fs3="fs pattern 4"
fs4="fs pattern 2"

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
            -color|--color) isColor=1;shift;;
            -ylabel|--ylabel) ylabel="$2";shift;;
            -title|--title) title="$2";shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -plot1|--plot1) isPlot1=1;;
            -multiplot|--multiplot) isMultiplot=1;;
            -showcase|--showcase) topt=$2; 
                if [ "${topt:0:1}" == "y" -o  "${topt:0:1}" == "Y" ]; then
                    isShowNumCase=1
                else
                    isShowNumCase=0
                fi
                shift
                ;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        dataFile=$1
    fi
    shift
done


if [ ! -f "$dataFile" ]; then 
    echo "Error! dataFile = \"$dataFile\" does not exist. Exit..."
    echo "$cmd"
    exit
fi

if [ "$outpath" == "" ]; then
    outpath=`dirname $dataFile`
elif [ ! -d "$outpath" ]; then
    mkdir -p $outpath
fi


stemname=${dataFile%.*}
dataFile1=$stemname.1.txt
dataFile2=$stemname.2.txt
dataFile3=$stemname.3.txt
dataFile4=$stemname.4.txt

if [ $isPlot1 -eq 1 ];then
    CreateDataFile1 "$dataFile" "$dataFile1"
    CreateDataFile2 "$dataFile" "$dataFile2"
    CreateDataFile3 "$dataFile" "$dataFile3"
    CreateDataFile4 "$dataFile" "$dataFile4"
fi


if [ $isShowNumCase -eq 1 ]; then
    MakePlot_shownumcase $dataFile
    if [ $isPlot1 -eq 1 ];then
        MakePlot_shownumcase $dataFile1
        MakePlot_shownumcase $dataFile2
    fi
else
    MakePlot_no $dataFile
    if [ $isPlot1 -eq 1 ];then
        for file in $dataFile1 $dataFile2 $dataFile3 $dataFile4; do
            MakePlot_no $file
            MakePlot_no_invert $file
            if [ $isMultiplot -eq 1 ];then
                MakePlot_no_invert_multiplot $file
            fi
        done
    fi
fi


